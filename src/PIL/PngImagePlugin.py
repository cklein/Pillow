#
# The Python Imaging Library.
# $Id$
#
# PNG support code
#
# See "PNG (Portable Network Graphics) Specification, version 1.0;
# W3C Recommendation", 1996-10-01, Thomas Boutell (ed.).
#
# history:
# 1996-05-06 fl   Created (couldn't resist it)
# 1996-12-14 fl   Upgraded, added read and verify support (0.2)
# 1996-12-15 fl   Separate PNG stream parser
# 1996-12-29 fl   Added write support, added getchunks
# 1996-12-30 fl   Eliminated circular references in decoder (0.3)
# 1998-07-12 fl   Read/write 16-bit images as mode I (0.4)
# 2001-02-08 fl   Added transparency support (from Zircon) (0.5)
# 2001-04-16 fl   Don't close data source in "open" method (0.6)
# 2004-02-24 fl   Don't even pretend to support interlaced files (0.7)
# 2004-08-31 fl   Do basic sanity check on chunk identifiers (0.8)
# 2004-09-20 fl   Added PngInfo chunk container
# 2004-12-18 fl   Added DPI read support (based on code by Niki Spahiev)
# 2008-08-13 fl   Added tRNS support for RGB images
# 2009-03-06 fl   Support for preserving ICC profiles (by Florian Hoech)
# 2009-03-08 fl   Added zTXT support (from Lowell Alleman)
# 2009-03-29 fl   Read interlaced PNG files (from Conrado Porto Lopes Gouvua)
#
# Copyright (c) 1997-2009 by Secret Labs AB
# Copyright (c) 1996 by Fredrik Lundh
#
# See the README file for information on usage and redistribution.
#

import itertools
import logging
import re
import struct
import warnings
import zlib

from . import Image, ImageChops, ImageFile, ImagePalette, ImageSequence
from ._binary import i8
from ._binary import i16be as i16
from ._binary import i32be as i32
from ._binary import o8
from ._binary import o16be as o16
from ._binary import o32be as o32

logger = logging.getLogger(__name__)

is_cid = re.compile(br"\w\w\w\w").match


_MAGIC = b"\211PNG\r\n\032\n"


_MODES = {
    # supported bits/color combinations, and corresponding modes/rawmodes
    # Greyscale
    (1, 0): ("1", "1"),
    (2, 0): ("L", "L;2"),
    (4, 0): ("L", "L;4"),
    (8, 0): ("L", "L"),
    (16, 0): ("I", "I;16B"),
    # Truecolour
    (8, 2): ("RGB", "RGB"),
    (16, 2): ("RGB", "RGB;16B"),
    # Indexed-colour
    (1, 3): ("P", "P;1"),
    (2, 3): ("P", "P;2"),
    (4, 3): ("P", "P;4"),
    (8, 3): ("P", "P"),
    # Greyscale with alpha
    (8, 4): ("LA", "LA"),
    (16, 4): ("RGBA", "LA;16B"),  # LA;16B->LA not yet available
    # Truecolour with alpha
    (8, 6): ("RGBA", "RGBA"),
    (16, 6): ("RGBA", "RGBA;16B"),
}


_simple_palette = re.compile(b"^\xff*\x00\xff*$")

MAX_TEXT_CHUNK = ImageFile.SAFEBLOCK
"""
Maximum decompressed size for a iTXt or zTXt chunk.
Eliminates decompression bombs where compressed chunks can expand 1000x.
See :ref:`Text in PNG File Format<png-text>`.
"""
MAX_TEXT_MEMORY = 64 * MAX_TEXT_CHUNK
"""
Set the maximum total text chunk size.
See :ref:`Text in PNG File Format<png-text>`.
"""


# APNG frame disposal modes
APNG_DISPOSE_OP_NONE = 0
"""
No disposal is done on this frame before rendering the next frame.
See :ref:`Saving APNG sequences<apng-saving>`.
"""
APNG_DISPOSE_OP_BACKGROUND = 1
"""
This frame’s modified region is cleared to fully transparent black before rendering
the next frame.
See :ref:`Saving APNG sequences<apng-saving>`.
"""
APNG_DISPOSE_OP_PREVIOUS = 2
"""
This frame’s modified region is reverted to the previous frame’s contents before
rendering the next frame.
See :ref:`Saving APNG sequences<apng-saving>`.
"""

# APNG frame blend modes
APNG_BLEND_OP_SOURCE = 0
"""
All color components of this frame, including alpha, overwrite the previous output
image contents.
See :ref:`Saving APNG sequences<apng-saving>`.
"""
APNG_BLEND_OP_OVER = 1
"""
This frame should be alpha composited with the previous output image contents.
See :ref:`Saving APNG sequences<apng-saving>`.
"""


def _safe_zlib_decompress(s):
    dobj = zlib.decompressobj()
    plaintext = dobj.decompress(s, MAX_TEXT_CHUNK)
    if dobj.unconsumed_tail:
        raise ValueError("Decompressed Data Too Large")
    return plaintext


def _crc32(data, seed=0):
    return zlib.crc32(data, seed) & 0xFFFFFFFF


# --------------------------------------------------------------------
# Support classes.  Suitable for PNG and related formats like MNG etc.


class ChunkStream:
    def __init__(self, fp):

        self.fp = fp
        self.queue = []

    def read(self):
        """Fetch a new chunk. Returns header information."""
        cid = None

        if self.queue:
            cid, pos, length = self.queue.pop()
            self.fp.seek(pos)
        else:
            s = self.fp.read(8)
            cid = s[4:]
            pos = self.fp.tell()
            length = i32(s)

        if not is_cid(cid):
            if not ImageFile.LOAD_TRUNCATED_IMAGES:
                raise SyntaxError("broken PNG file (chunk %s)" % repr(cid))

        return cid, pos, length

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        self.queue = self.crc = self.fp = None

    def push(self, cid, pos, length):

        self.queue.append((cid, pos, length))

    def call(self, cid, pos, length):
        """Call the appropriate chunk handler"""

        logger.debug("STREAM %r %s %s", cid, pos, length)
        return getattr(self, "chunk_" + cid.decode("ascii"))(pos, length)

    def crc(self, cid, data):
        """Read and verify checksum"""

        # Skip CRC checks for ancillary chunks if allowed to load truncated
        # images
        # 5th byte of first char is 1 [specs, section 5.4]
        if ImageFile.LOAD_TRUNCATED_IMAGES and (i8(cid[0]) >> 5 & 1):
            self.crc_skip(cid, data)
            return

        try:
            crc1 = _crc32(data, _crc32(cid))
            crc2 = i32(self.fp.read(4))
            if crc1 != crc2:
                raise SyntaxError("broken PNG file (bad header checksum in %r)" % cid)
        except struct.error as e:
            raise SyntaxError(
                "broken PNG file (incomplete checksum in %r)" % cid
            ) from e

    def crc_skip(self, cid, data):
        """Read checksum.  Used if the C module is not present"""

        self.fp.read(4)

    def verify(self, endchunk=b"IEND"):

        # Simple approach; just calculate checksum for all remaining
        # blocks.  Must be called directly after open.

        cids = []

        while True:
            try:
                cid, pos, length = self.read()
            except struct.error as e:
                raise OSError("truncated PNG file") from e

            if cid == endchunk:
                break
            self.crc(cid, ImageFile._safe_read(self.fp, length))
            cids.append(cid)

        return cids


class iTXt(str):
    """
    Subclass of string to allow iTXt chunks to look like strings while
    keeping their extra information

    """

    @staticmethod
    def __new__(cls, text, lang=None, tkey=None):
        """
        :param cls: the class to use when creating the instance
        :param text: value for this key
        :param lang: language code
        :param tkey: UTF-8 version of the key name
        """

        self = str.__new__(cls, text)
        self.lang = lang
        self.tkey = tkey
        return self


class PngInfo:
    """
    PNG chunk container (for use with save(pnginfo=))

    """

    def __init__(self):
        self.chunks = []

    def add(self, cid, data):
        """Appends an arbitrary chunk. Use with caution.

        :param cid: a byte string, 4 bytes long.
        :param data: a byte string of the encoded data

        """

        self.chunks.append((cid, data))

    def add_itxt(self, key, value, lang="", tkey="", zip=False):
        """Appends an iTXt chunk.

        :param key: latin-1 encodable text key name
        :param value: value for this key
        :param lang: language code
        :param tkey: UTF-8 version of the key name
        :param zip: compression flag

        """

        if not isinstance(key, bytes):
            key = key.encode("latin-1", "strict")
        if not isinstance(value, bytes):
            value = value.encode("utf-8", "strict")
        if not isinstance(lang, bytes):
            lang = lang.encode("utf-8", "strict")
        if not isinstance(tkey, bytes):
            tkey = tkey.encode("utf-8", "strict")

        if zip:
            self.add(
                b"iTXt",
                key + b"\0\x01\0" + lang + b"\0" + tkey + b"\0" + zlib.compress(value),
            )
        else:
            self.add(b"iTXt", key + b"\0\0\0" + lang + b"\0" + tkey + b"\0" + value)

    def add_text(self, key, value, zip=False):
        """Appends a text chunk.

        :param key: latin-1 encodable text key name
        :param value: value for this key, text or an
           :py:class:`PIL.PngImagePlugin.iTXt` instance
        :param zip: compression flag

        """
        if isinstance(value, iTXt):
            return self.add_itxt(key, value, value.lang, value.tkey, zip=zip)

        # The tEXt chunk stores latin-1 text
        if not isinstance(value, bytes):
            try:
                value = value.encode("latin-1", "strict")
            except UnicodeError:
                return self.add_itxt(key, value, zip=zip)

        if not isinstance(key, bytes):
            key = key.encode("latin-1", "strict")

        if zip:
            self.add(b"zTXt", key + b"\0\0" + zlib.compress(value))
        else:
            self.add(b"tEXt", key + b"\0" + value)


# --------------------------------------------------------------------
# PNG image stream (IHDR/IEND)


class PngStream(ChunkStream):
    def __init__(self, fp):
        super().__init__(fp)

        # local copies of Image attributes
        self.im_info = {}
        self.im_text = {}
        self.im_size = (0, 0)
        self.im_mode = None
        self.im_tile = None
        self.im_palette = None
        self.im_custom_mimetype = None
        self.im_n_frames = None
        self._seq_num = None
        self.rewind_state = None

        self.text_memory = 0

    def check_text_memory(self, chunklen):
        self.text_memory += chunklen
        if self.text_memory > MAX_TEXT_MEMORY:
            raise ValueError(
                "Too much memory used in text chunks: %s>MAX_TEXT_MEMORY"
                % self.text_memory
            )

    def save_rewind(self):
        self.rewind_state = {
            "info": self.im_info.copy(),
            "tile": self.im_tile,
            "seq_num": self._seq_num,
        }

    def rewind(self):
        self.im_info = self.rewind_state["info"]
        self.im_tile = self.rewind_state["tile"]
        self._seq_num = self.rewind_state["seq_num"]

    def chunk_iCCP(self, pos, length):

        # ICC profile
        s = ImageFile._safe_read(self.fp, length)
        # according to PNG spec, the iCCP chunk contains:
        # Profile name  1-79 bytes (character string)
        # Null separator        1 byte (null character)
        # Compression method    1 byte (0)
        # Compressed profile    n bytes (zlib with deflate compression)
        i = s.find(b"\0")
        logger.debug("iCCP profile name %r", s[:i])
        logger.debug("Compression method %s", i8(s[i]))
        comp_method = i8(s[i])
        if comp_method != 0:
            raise SyntaxError(
                "Unknown compression method %s in iCCP chunk" % comp_method
            )
        try:
            icc_profile = _safe_zlib_decompress(s[i + 2 :])
        except ValueError:
            if ImageFile.LOAD_TRUNCATED_IMAGES:
                icc_profile = None
            else:
                raise
        except zlib.error:
            icc_profile = None  # FIXME
        self.im_info["icc_profile"] = icc_profile
        return s

    def chunk_IHDR(self, pos, length):

        # image header
        s = ImageFile._safe_read(self.fp, length)
        self.im_size = i32(s), i32(s[4:])
        try:
            self.im_mode, self.im_rawmode = _MODES[(i8(s[8]), i8(s[9]))]
        except Exception:
            pass
        if i8(s[12]):
            self.im_info["interlace"] = 1
        if i8(s[11]):
            raise SyntaxError("unknown filter category")
        return s

    def chunk_IDAT(self, pos, length):

        # image data
        if "bbox" in self.im_info:
            tile = [("zip", self.im_info["bbox"], pos, self.im_rawmode)]
        else:
            if self.im_n_frames is not None:
                self.im_info["default_image"] = True
            tile = [("zip", (0, 0) + self.im_size, pos, self.im_rawmode)]
        self.im_tile = tile
        self.im_idat = length
        raise EOFError

    def chunk_IEND(self, pos, length):

        # end of PNG image
        raise EOFError

    def chunk_PLTE(self, pos, length):

        # palette
        s = ImageFile._safe_read(self.fp, length)
        if self.im_mode == "P":
            self.im_palette = "RGB", s
        return s

    def chunk_tRNS(self, pos, length):

        # transparency
        s = ImageFile._safe_read(self.fp, length)
        if self.im_mode == "P":
            if _simple_palette.match(s):
                # tRNS contains only one full-transparent entry,
                # other entries are full opaque
                i = s.find(b"\0")
                if i >= 0:
                    self.im_info["transparency"] = i
            else:
                # otherwise, we have a byte string with one alpha value
                # for each palette entry
                self.im_info["transparency"] = s
        elif self.im_mode in ("1", "L", "I"):
            self.im_info["transparency"] = i16(s)
        elif self.im_mode == "RGB":
            self.im_info["transparency"] = i16(s), i16(s[2:]), i16(s[4:])
        return s

    def chunk_gAMA(self, pos, length):
        # gamma setting
        s = ImageFile._safe_read(self.fp, length)
        self.im_info["gamma"] = i32(s) / 100000.0
        return s

    def chunk_cHRM(self, pos, length):
        # chromaticity, 8 unsigned ints, actual value is scaled by 100,000
        # WP x,y, Red x,y, Green x,y Blue x,y

        s = ImageFile._safe_read(self.fp, length)
        raw_vals = struct.unpack(">%dI" % (len(s) // 4), s)
        self.im_info["chromaticity"] = tuple(elt / 100000.0 for elt in raw_vals)
        return s

    def chunk_sRGB(self, pos, length):
        # srgb rendering intent, 1 byte
        # 0 perceptual
        # 1 relative colorimetric
        # 2 saturation
        # 3 absolute colorimetric

        s = ImageFile._safe_read(self.fp, length)
        self.im_info["srgb"] = i8(s)
        return s

    def chunk_pHYs(self, pos, length):

        # pixels per unit
        s = ImageFile._safe_read(self.fp, length)
        px, py = i32(s), i32(s[4:])
        unit = i8(s[8])
        if unit == 1:  # meter
            dpi = int(px * 0.0254 + 0.5), int(py * 0.0254 + 0.5)
            self.im_info["dpi"] = dpi
        elif unit == 0:
            self.im_info["aspect"] = px, py
        return s

    def chunk_tEXt(self, pos, length):

        # text
        s = ImageFile._safe_read(self.fp, length)
        try:
            k, v = s.split(b"\0", 1)
        except ValueError:
            # fallback for broken tEXt tags
            k = s
            v = b""
        if k:
            k = k.decode("latin-1", "strict")
            v_str = v.decode("latin-1", "replace")

            self.im_info[k] = v if k == "exif" else v_str
            self.im_text[k] = v_str
            self.check_text_memory(len(v_str))

        return s

    def chunk_zTXt(self, pos, length):

        # compressed text
        s = ImageFile._safe_read(self.fp, length)
        try:
            k, v = s.split(b"\0", 1)
        except ValueError:
            k = s
            v = b""
        if v:
            comp_method = i8(v[0])
        else:
            comp_method = 0
        if comp_method != 0:
            raise SyntaxError(
                "Unknown compression method %s in zTXt chunk" % comp_method
            )
        try:
            v = _safe_zlib_decompress(v[1:])
        except ValueError:
            if ImageFile.LOAD_TRUNCATED_IMAGES:
                v = b""
            else:
                raise
        except zlib.error:
            v = b""

        if k:
            k = k.decode("latin-1", "strict")
            v = v.decode("latin-1", "replace")

            self.im_info[k] = self.im_text[k] = v
            self.check_text_memory(len(v))

        return s

    def chunk_iTXt(self, pos, length):

        # international text
        r = s = ImageFile._safe_read(self.fp, length)
        try:
            k, r = r.split(b"\0", 1)
        except ValueError:
            return s
        if len(r) < 2:
            return s
        cf, cm, r = i8(r[0]), i8(r[1]), r[2:]
        try:
            lang, tk, v = r.split(b"\0", 2)
        except ValueError:
            return s
        if cf != 0:
            if cm == 0:
                try:
                    v = _safe_zlib_decompress(v)
                except ValueError:
                    if ImageFile.LOAD_TRUNCATED_IMAGES:
                        return s
                    else:
                        raise
                except zlib.error:
                    return s
            else:
                return s
        try:
            k = k.decode("latin-1", "strict")
            lang = lang.decode("utf-8", "strict")
            tk = tk.decode("utf-8", "strict")
            v = v.decode("utf-8", "strict")
        except UnicodeError:
            return s

        self.im_info[k] = self.im_text[k] = iTXt(v, lang, tk)
        self.check_text_memory(len(v))

        return s

    def chunk_eXIf(self, pos, length):
        s = ImageFile._safe_read(self.fp, length)
        self.im_info["exif"] = b"Exif\x00\x00" + s
        return s

    # APNG chunks
    def chunk_acTL(self, pos, length):
        s = ImageFile._safe_read(self.fp, length)
        if self.im_n_frames is not None:
            self.im_n_frames = None
            warnings.warn("Invalid APNG, will use default PNG image if possible")
            return s
        n_frames = i32(s)
        if n_frames == 0 or n_frames > 0x80000000:
            warnings.warn("Invalid APNG, will use default PNG image if possible")
            return s
        self.im_n_frames = n_frames
        self.im_info["loop"] = i32(s[4:])
        self.im_custom_mimetype = "image/apng"
        return s

    def chunk_fcTL(self, pos, length):
        s = ImageFile._safe_read(self.fp, length)
        seq = i32(s)
        if (self._seq_num is None and seq != 0) or (
            self._seq_num is not None and self._seq_num != seq - 1
        ):
            raise SyntaxError("APNG contains frame sequence errors")
        self._seq_num = seq
        width, height = i32(s[4:]), i32(s[8:])
        px, py = i32(s[12:]), i32(s[16:])
        im_w, im_h = self.im_size
        if px + width > im_w or py + height > im_h:
            raise SyntaxError("APNG contains invalid frames")
        self.im_info["bbox"] = (px, py, px + width, py + height)
        delay_num, delay_den = i16(s[20:]), i16(s[22:])
        if delay_den == 0:
            delay_den = 100
        self.im_info["duration"] = float(delay_num) / float(delay_den) * 1000
        self.im_info["disposal"] = i8(s[24])
        self.im_info["blend"] = i8(s[25])
        return s

    def chunk_fdAT(self, pos, length):
        s = ImageFile._safe_read(self.fp, 4)
        seq = i32(s)
        if self._seq_num != seq - 1:
            raise SyntaxError("APNG contains frame sequence errors")
        self._seq_num = seq
        return self.chunk_IDAT(pos + 4, length - 4)


# --------------------------------------------------------------------
# PNG reader


def _accept(prefix):
    return prefix[:8] == _MAGIC


##
# Image plugin for PNG images.


class PngImageFile(ImageFile.ImageFile):

    format = "PNG"
    format_description = "Portable network graphics"

    def _open(self):

        if not _accept(self.fp.read(8)):
            raise SyntaxError("not a PNG file")
        self.__fp = self.fp
        self.__frame = 0

        #
        # Parse headers up to the first IDAT or fDAT chunk

        self.png = PngStream(self.fp)

        while True:

            #
            # get next chunk

            cid, pos, length = self.png.read()

            try:
                s = self.png.call(cid, pos, length)
            except EOFError:
                break
            except AttributeError:
                logger.debug("%r %s %s (unknown)", cid, pos, length)
                s = ImageFile._safe_read(self.fp, length)

            self.png.crc(cid, s)

        #
        # Copy relevant attributes from the PngStream.  An alternative
        # would be to let the PngStream class modify these attributes
        # directly, but that introduces circular references which are
        # difficult to break if things go wrong in the decoder...
        # (believe me, I've tried ;-)

        self.mode = self.png.im_mode
        self._size = self.png.im_size
        self.info = self.png.im_info
        self._text = None
        self.tile = self.png.im_tile
        self.custom_mimetype = self.png.im_custom_mimetype
        self.n_frames = self.png.im_n_frames or 1
        self.default_image = self.info.get("default_image", False)

        if self.png.im_palette:
            rawmode, data = self.png.im_palette
            self.palette = ImagePalette.raw(rawmode, data)

        if cid == b"fdAT":
            self.__prepare_idat = length - 4
        else:
            self.__prepare_idat = length  # used by load_prepare()

        if self.png.im_n_frames is not None:
            self._close_exclusive_fp_after_loading = False
            self.png.save_rewind()
            self.__rewind_idat = self.__prepare_idat
            self.__rewind = self.__fp.tell()
            if self.default_image:
                # IDAT chunk contains default image and not first animation frame
                self.n_frames += 1
            self._seek(0)
        self.is_animated = self.n_frames > 1

    @property
    def text(self):
        # experimental
        if self._text is None:
            # iTxt, tEXt and zTXt chunks may appear at the end of the file
            # So load the file to ensure that they are read
            if self.is_animated:
                frame = self.__frame
                # for APNG, seek to the final frame before loading
                self.seek(self.n_frames - 1)
            self.load()
            if self.is_animated:
                self.seek(frame)
        return self._text

    def verify(self):
        """Verify PNG file"""

        if self.fp is None:
            raise RuntimeError("verify must be called directly after open")

        # back up to beginning of IDAT block
        self.fp.seek(self.tile[0][2] - 8)

        self.png.verify()
        self.png.close()

        if self._exclusive_fp:
            self.fp.close()
        self.fp = None

    def seek(self, frame):
        if not self._seek_check(frame):
            return
        if frame < self.__frame:
            self._seek(0, True)

        last_frame = self.__frame
        for f in range(self.__frame + 1, frame + 1):
            try:
                self._seek(f)
            except EOFError as e:
                self.seek(last_frame)
                raise EOFError("no more images in APNG file") from e

    def _seek(self, frame, rewind=False):
        if frame == 0:
            if rewind:
                self.__fp.seek(self.__rewind)
                self.png.rewind()
                self.__prepare_idat = self.__rewind_idat
                self.im = None
                if self.pyaccess:
                    self.pyaccess = None
                self.info = self.png.im_info
                self.tile = self.png.im_tile
                self.fp = self.__fp
            self._prev_im = None
            self.dispose = None
            self.default_image = self.info.get("default_image", False)
            self.dispose_op = self.info.get("disposal")
            self.blend_op = self.info.get("blend")
            self.dispose_extent = self.info.get("bbox")
            self.__frame = 0
            return
        else:
            if frame != self.__frame + 1:
                raise ValueError("cannot seek to frame %d" % frame)

        # ensure previous frame was loaded
        self.load()

        self.fp = self.__fp

        # advance to the next frame
        if self.__prepare_idat:
            ImageFile._safe_read(self.fp, self.__prepare_idat)
            self.__prepare_idat = 0
        frame_start = False
        while True:
            self.fp.read(4)  # CRC

            try:
                cid, pos, length = self.png.read()
            except (struct.error, SyntaxError):
                break

            if cid == b"IEND":
                raise EOFError("No more images in APNG file")
            if cid == b"fcTL":
                if frame_start:
                    # there must be at least one fdAT chunk between fcTL chunks
                    raise SyntaxError("APNG missing frame data")
                frame_start = True

            try:
                self.png.call(cid, pos, length)
            except UnicodeDecodeError:
                break
            except EOFError:
                if cid == b"fdAT":
                    length -= 4
                    if frame_start:
                        self.__prepare_idat = length
                        break
                ImageFile._safe_read(self.fp, length)
            except AttributeError:
                logger.debug("%r %s %s (unknown)", cid, pos, length)
                ImageFile._safe_read(self.fp, length)

        self.__frame = frame
        self.tile = self.png.im_tile
        self.dispose_op = self.info.get("disposal")
        self.blend_op = self.info.get("blend")
        self.dispose_extent = self.info.get("bbox")

        if not self.tile:
            raise EOFError

    def tell(self):
        return self.__frame

    def load_prepare(self):
        """internal: prepare to read PNG file"""

        if self.info.get("interlace"):
            self.decoderconfig = self.decoderconfig + (1,)

        self.__idat = self.__prepare_idat  # used by load_read()
        ImageFile.ImageFile.load_prepare(self)

    def load_read(self, read_bytes):
        """internal: read more image data"""

        while self.__idat == 0:
            # end of chunk, skip forward to next one

            self.fp.read(4)  # CRC

            cid, pos, length = self.png.read()

            if cid not in [b"IDAT", b"DDAT", b"fdAT"]:
                self.png.push(cid, pos, length)
                return b""

            if cid == b"fdAT":
                try:
                    self.png.call(cid, pos, length)
                except EOFError:
                    pass
                self.__idat = length - 4  # sequence_num has already been read
            else:
                self.__idat = length  # empty chunks are allowed

        # read more data from this chunk
        if read_bytes <= 0:
            read_bytes = self.__idat
        else:
            read_bytes = min(read_bytes, self.__idat)

        self.__idat = self.__idat - read_bytes

        return self.fp.read(read_bytes)

    def load_end(self):
        """internal: finished reading image data"""
        while True:
            self.fp.read(4)  # CRC

            try:
                cid, pos, length = self.png.read()
            except (struct.error, SyntaxError):
                break

            if cid == b"IEND":
                break
            elif cid == b"fcTL" and self.is_animated:
                # start of the next frame, stop reading
                self.__prepare_idat = 0
                self.png.push(cid, pos, length)
                break

            try:
                self.png.call(cid, pos, length)
            except UnicodeDecodeError:
                break
            except EOFError:
                if cid == b"fdAT":
                    length -= 4
                ImageFile._safe_read(self.fp, length)
            except AttributeError:
                logger.debug("%r %s %s (unknown)", cid, pos, length)
                ImageFile._safe_read(self.fp, length)
        self._text = self.png.im_text
        if not self.is_animated:
            self.png.close()
            self.png = None
        else:
            # setup frame disposal (actual disposal done when needed in _seek())
            if self._prev_im is None and self.dispose_op == APNG_DISPOSE_OP_PREVIOUS:
                self.dispose_op = APNG_DISPOSE_OP_BACKGROUND

            if self.dispose_op == APNG_DISPOSE_OP_PREVIOUS:
                dispose = self._prev_im.copy()
                dispose = self._crop(dispose, self.dispose_extent)
            elif self.dispose_op == APNG_DISPOSE_OP_BACKGROUND:
                dispose = Image.core.fill(self.im.mode, self.size)
                dispose = self._crop(dispose, self.dispose_extent)
            else:
                dispose = None

            if self._prev_im and self.blend_op == APNG_BLEND_OP_OVER:
                updated = self._crop(self.im, self.dispose_extent)
                self._prev_im.paste(
                    updated, self.dispose_extent, updated.convert("RGBA")
                )
                self.im = self._prev_im
                if self.pyaccess:
                    self.pyaccess = None
            self._prev_im = self.im.copy()

            if dispose:
                self._prev_im.paste(dispose, self.dispose_extent)

    def _getexif(self):
        if "exif" not in self.info:
            self.load()
        if "exif" not in self.info and "Raw profile type exif" not in self.info:
            return None
        return dict(self.getexif())

    def getexif(self):
        if "exif" not in self.info:
            self.load()

        return super().getexif()

    def _close__fp(self):
        try:
            if self.__fp != self.fp:
                self.__fp.close()
        except AttributeError:
            pass
        finally:
            self.__fp = None


# --------------------------------------------------------------------
# PNG writer

_OUTMODES = {
    # supported PIL modes, and corresponding rawmodes/bits/color combinations
    "1": ("1", b"\x01\x00"),
    "L;1": ("L;1", b"\x01\x00"),
    "L;2": ("L;2", b"\x02\x00"),
    "L;4": ("L;4", b"\x04\x00"),
    "L": ("L", b"\x08\x00"),
    "LA": ("LA", b"\x08\x04"),
    "I": ("I;16B", b"\x10\x00"),
    "I;16": ("I;16B", b"\x10\x00"),
    "P;1": ("P;1", b"\x01\x03"),
    "P;2": ("P;2", b"\x02\x03"),
    "P;4": ("P;4", b"\x04\x03"),
    "P": ("P", b"\x08\x03"),
    "RGB": ("RGB", b"\x08\x02"),
    "RGBA": ("RGBA", b"\x08\x06"),
}


def putchunk(fp, cid, *data):
    """Write a PNG chunk (including CRC field)"""

    data = b"".join(data)

    fp.write(o32(len(data)) + cid)
    fp.write(data)
    crc = _crc32(data, _crc32(cid))
    fp.write(o32(crc))


class _idat:
    # wrap output from the encoder in IDAT chunks

    def __init__(self, fp, chunk):
        self.fp = fp
        self.chunk = chunk

    def write(self, data):
        self.chunk(self.fp, b"IDAT", data)


class _fdat:
    # wrap encoder output in fdAT chunks

    def __init__(self, fp, chunk, seq_num):
        self.fp = fp
        self.chunk = chunk
        self.seq_num = seq_num

    def write(self, data):
        self.chunk(self.fp, b"fdAT", o32(self.seq_num), data)
        self.seq_num += 1


def _write_multiple_frames(im, fp, chunk, rawmode):
    default_image = im.encoderinfo.get("default_image", im.info.get("default_image"))
    duration = im.encoderinfo.get("duration", im.info.get("duration", 0))
    loop = im.encoderinfo.get("loop", im.info.get("loop", 0))
    disposal = im.encoderinfo.get("disposal", im.info.get("disposal"))
    blend = im.encoderinfo.get("blend", im.info.get("blend"))

    if default_image:
        chain = itertools.chain(im.encoderinfo.get("append_images", []))
    else:
        chain = itertools.chain([im], im.encoderinfo.get("append_images", []))

    im_frames = []
    frame_count = 0
    for im_seq in chain:
        for im_frame in ImageSequence.Iterator(im_seq):
            im_frame = im_frame.copy()
            if im_frame.mode != im.mode:
                if im.mode == "P":
                    im_frame = im_frame.convert(im.mode, palette=im.palette)
                else:
                    im_frame = im_frame.convert(im.mode)
            encoderinfo = im.encoderinfo.copy()
            if isinstance(duration, (list, tuple)):
                encoderinfo["duration"] = duration[frame_count]
            if isinstance(disposal, (list, tuple)):
                encoderinfo["disposal"] = disposal[frame_count]
            if isinstance(blend, (list, tuple)):
                encoderinfo["blend"] = blend[frame_count]
            frame_count += 1

            if im_frames:
                previous = im_frames[-1]
                prev_disposal = previous["encoderinfo"].get("disposal")
                prev_blend = previous["encoderinfo"].get("blend")
                if prev_disposal == APNG_DISPOSE_OP_PREVIOUS and len(im_frames) < 2:
                    prev_disposal = APNG_DISPOSE_OP_BACKGROUND

                if prev_disposal == APNG_DISPOSE_OP_BACKGROUND:
                    base_im = previous["im"]
                    dispose = Image.core.fill("RGBA", im.size, (0, 0, 0, 0))
                    bbox = previous["bbox"]
                    if bbox:
                        dispose = dispose.crop(bbox)
                    else:
                        bbox = (0, 0) + im.size
                    base_im.paste(dispose, bbox)
                elif prev_disposal == APNG_DISPOSE_OP_PREVIOUS:
                    base_im = im_frames[-2]["im"]
                else:
                    base_im = previous["im"]
                delta = ImageChops.subtract_modulo(
                    im_frame.convert("RGB"), base_im.convert("RGB")
                )
                bbox = delta.getbbox()
                if (
                    not bbox
                    and prev_disposal == encoderinfo.get("disposal")
                    and prev_blend == encoderinfo.get("blend")
                ):
                    duration = encoderinfo.get("duration", 0)
                    if duration:
                        if "duration" in previous["encoderinfo"]:
                            previous["encoderinfo"]["duration"] += duration
                        else:
                            previous["encoderinfo"]["duration"] = duration
                    continue
            else:
                bbox = None
            im_frames.append({"im": im_frame, "bbox": bbox, "encoderinfo": encoderinfo})

    # animation control
    chunk(
        fp, b"acTL", o32(len(im_frames)), o32(loop),  # 0: num_frames  # 4: num_plays
    )

    # default image IDAT (if it exists)
    if default_image:
        ImageFile._save(im, _idat(fp, chunk), [("zip", (0, 0) + im.size, 0, rawmode)])

    seq_num = 0
    for frame, frame_data in enumerate(im_frames):
        im_frame = frame_data["im"]
        if not frame_data["bbox"]:
            bbox = (0, 0) + im_frame.size
        else:
            bbox = frame_data["bbox"]
            im_frame = im_frame.crop(bbox)
        size = im_frame.size
        duration = int(round(frame_data["encoderinfo"].get("duration", 0)))
        disposal = frame_data["encoderinfo"].get("disposal", APNG_DISPOSE_OP_NONE)
        blend = frame_data["encoderinfo"].get("blend", APNG_BLEND_OP_SOURCE)
        # frame control
        chunk(
            fp,
            b"fcTL",
            o32(seq_num),  # sequence_number
            o32(size[0]),  # width
            o32(size[1]),  # height
            o32(bbox[0]),  # x_offset
            o32(bbox[1]),  # y_offset
            o16(duration),  # delay_numerator
            o16(1000),  # delay_denominator
            o8(disposal),  # dispose_op
            o8(blend),  # blend_op
        )
        seq_num += 1
        # frame data
        if frame == 0 and not default_image:
            # first frame must be in IDAT chunks for backwards compatibility
            ImageFile._save(
                im_frame,
                _idat(fp, chunk),
                [("zip", (0, 0) + im_frame.size, 0, rawmode)],
            )
        else:
            fdat_chunks = _fdat(fp, chunk, seq_num)
            ImageFile._save(
                im_frame, fdat_chunks, [("zip", (0, 0) + im_frame.size, 0, rawmode)],
            )
            seq_num = fdat_chunks.seq_num


def _save_all(im, fp, filename):
    _save(im, fp, filename, save_all=True)


def _save(im, fp, filename, chunk=putchunk, save_all=False):
    # save an image to disk (called by the save method)

    mode = im.mode

    if mode == "P":

        #
        # attempt to minimize storage requirements for palette images
        if "bits" in im.encoderinfo:
            # number of bits specified by user
            colors = 1 << im.encoderinfo["bits"]
        else:
            # check palette contents
            if im.palette:
                colors = max(min(len(im.palette.getdata()[1]) // 3, 256), 2)
            else:
                colors = 256

        if colors <= 2:
            bits = 1
        elif colors <= 4:
            bits = 2
        elif colors <= 16:
            bits = 4
        else:
            bits = 8
        if bits != 8:
            mode = "%s;%d" % (mode, bits)

    # encoder options
    im.encoderconfig = (
        im.encoderinfo.get("optimize", False),
        im.encoderinfo.get("compress_level", -1),
        im.encoderinfo.get("compress_type", -1),
        im.encoderinfo.get("dictionary", b""),
    )

    # get the corresponding PNG mode
    try:
        rawmode, mode = _OUTMODES[mode]
    except KeyError as e:
        raise OSError("cannot write mode %s as PNG" % mode) from e

    #
    # write minimal PNG file

    fp.write(_MAGIC)

    chunk(
        fp,
        b"IHDR",
        o32(im.size[0]),  # 0: size
        o32(im.size[1]),
        mode,  # 8: depth/type
        b"\0",  # 10: compression
        b"\0",  # 11: filter category
        b"\0",  # 12: interlace flag
    )

    chunks = [b"cHRM", b"gAMA", b"sBIT", b"sRGB", b"tIME"]

    icc = im.encoderinfo.get("icc_profile", im.info.get("icc_profile"))
    if icc:
        # ICC profile
        # according to PNG spec, the iCCP chunk contains:
        # Profile name  1-79 bytes (character string)
        # Null separator        1 byte (null character)
        # Compression method    1 byte (0)
        # Compressed profile    n bytes (zlib with deflate compression)
        name = b"ICC Profile"
        data = name + b"\0\0" + zlib.compress(icc)
        chunk(fp, b"iCCP", data)

        # You must either have sRGB or iCCP.
        # Disallow sRGB chunks when an iCCP-chunk has been emitted.
        chunks.remove(b"sRGB")

    info = im.encoderinfo.get("pnginfo")
    if info:
        chunks_multiple_allowed = [b"sPLT", b"iTXt", b"tEXt", b"zTXt"]
        for cid, data in info.chunks:
            if cid in chunks:
                chunks.remove(cid)
                chunk(fp, cid, data)
            elif cid in chunks_multiple_allowed:
                chunk(fp, cid, data)

    if im.mode == "P":
        palette_byte_number = (2 ** bits) * 3
        palette_bytes = im.im.getpalette("RGB")[:palette_byte_number]
        while len(palette_bytes) < palette_byte_number:
            palette_bytes += b"\0"
        chunk(fp, b"PLTE", palette_bytes)

    transparency = im.encoderinfo.get("transparency", im.info.get("transparency", None))

    if transparency or transparency == 0:
        if im.mode == "P":
            # limit to actual palette size
            alpha_bytes = 2 ** bits
            if isinstance(transparency, bytes):
                chunk(fp, b"tRNS", transparency[:alpha_bytes])
            else:
                transparency = max(0, min(255, transparency))
                alpha = b"\xFF" * transparency + b"\0"
                chunk(fp, b"tRNS", alpha[:alpha_bytes])
        elif im.mode in ("1", "L", "I"):
            transparency = max(0, min(65535, transparency))
            chunk(fp, b"tRNS", o16(transparency))
        elif im.mode == "RGB":
            red, green, blue = transparency
            chunk(fp, b"tRNS", o16(red) + o16(green) + o16(blue))
        else:
            if "transparency" in im.encoderinfo:
                # don't bother with transparency if it's an RGBA
                # and it's in the info dict. It's probably just stale.
                raise OSError("cannot use transparency for this mode")
    else:
        if im.mode == "P" and im.im.getpalettemode() == "RGBA":
            alpha = im.im.getpalette("RGBA", "A")
            alpha_bytes = 2 ** bits
            chunk(fp, b"tRNS", alpha[:alpha_bytes])

    dpi = im.encoderinfo.get("dpi")
    if dpi:
        chunk(
            fp,
            b"pHYs",
            o32(int(dpi[0] / 0.0254 + 0.5)),
            o32(int(dpi[1] / 0.0254 + 0.5)),
            b"\x01",
        )

    if info:
        chunks = [b"bKGD", b"hIST"]
        for cid, data in info.chunks:
            if cid in chunks:
                chunks.remove(cid)
                chunk(fp, cid, data)

    exif = im.encoderinfo.get("exif", im.info.get("exif"))
    if exif:
        if isinstance(exif, Image.Exif):
            exif = exif.tobytes(8)
        if exif.startswith(b"Exif\x00\x00"):
            exif = exif[6:]
        chunk(fp, b"eXIf", exif)

    if save_all:
        _write_multiple_frames(im, fp, chunk, rawmode)
    else:
        ImageFile._save(im, _idat(fp, chunk), [("zip", (0, 0) + im.size, 0, rawmode)])

    chunk(fp, b"IEND", b"")

    if hasattr(fp, "flush"):
        fp.flush()


# --------------------------------------------------------------------
# PNG chunk converter


def getchunks(im, **params):
    """Return a list of PNG chunks representing this image."""

    class collector:
        data = []

        def write(self, data):
            pass

        def append(self, chunk):
            self.data.append(chunk)

    def append(fp, cid, *data):
        data = b"".join(data)
        crc = o32(_crc32(data, _crc32(cid)))
        fp.append((cid, data, crc))

    fp = collector()

    try:
        im.encoderinfo = params
        _save(im, fp, None, append)
    finally:
        del im.encoderinfo

    return fp.data


# --------------------------------------------------------------------
# Registry

Image.register_open(PngImageFile.format, PngImageFile, _accept)
Image.register_save(PngImageFile.format, _save)
Image.register_save_all(PngImageFile.format, _save_all)

Image.register_extensions(PngImageFile.format, [".png", ".apng"])

Image.register_mime(PngImageFile.format, "image/png")
