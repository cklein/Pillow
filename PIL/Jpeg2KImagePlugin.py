#
# The Python Imaging Library
# $Id$
#
# JPEG2000 file handling
#
# History:
# 2014-03-12 ajh  Created
#
# Copyright (c) 2014 Coriolis Systems Limited
# Copyright (c) 2014 Alastair Houghton
#
# See the README file for information on usage and redistribution.
#

__version__ = "0.1"

from PIL import Image, ImageFile, _binary
import struct
import os
import io

def _parse_codestream(fp):
    """Parse the JPEG 2000 codestream to extract the size and component
    count from the SIZ marker segment, returning a PIL (size, mode) tuple."""
    
    hdr = fp.read(2)
    lsiz = struct.unpack('>H', hdr)[0]
    siz = hdr + fp.read(lsiz - 2)
    lsiz, rsiz, xsiz, ysiz, xosiz, yosiz, xtsiz, ytsiz, \
    xtosiz, ytosiz, csiz \
        = struct.unpack('>HHIIIIIIIIH', siz[:38])
    ssiz = [None]*csiz
    xrsiz = [None]*csiz
    yrsiz = [None]*csiz
    for i in range(csiz):
        ssiz[i], xrsiz[i], yrsiz[i] \
            = struct.unpack('>BBB', siz[36 + 3 * i:39 + 3 * i])

    size = (xsiz - xosiz, ysiz - yosiz)
    if csiz == 1:
        mode = 'L'
    elif csiz == 2:
        mode = 'LA'
    elif csiz == 3:
        mode = 'RGB'
    elif csiz == 4:
        mode == 'RGBA'
    else:
        mode = None
    
    return (size, mode)

def _parse_jp2_header(fp):
    """Parse the JP2 header box to extract size, component count and
    color space information, returning a PIL (size, mode) tuple."""
    
    # Find the JP2 header box
    header = None
    while True:
        lbox, tbox = struct.unpack('>I4s', fp.read(8))
        if lbox == 1:
            lbox = struct.unpack('>Q', fp.read(8))[0]
            hlen = 16
        else:
            hlen = 8

        if tbox == b'jp2h':
            header = fp.read(lbox - hlen)
            break
        else:
            fp.seek(lbox - hlen, os.SEEK_CUR)

    if header is None:
        raise SyntaxError('could not find JP2 header')

    size = None
    mode = None
    
    hio = io.BytesIO(header)
    while True:
        lbox, tbox = struct.unpack('>I4s', hio.read(8))
        if lbox == 1:
            lbox = struct.unpack('>Q', hio.read(8))[0]
            hlen = 16
        else:
            hlen = 8

        content = hio.read(lbox - hlen)

        if tbox == b'ihdr':
            height, width, nc, bpc, c, unkc, ipr \
              = struct.unpack('>IIHBBBB', content)
            size = (width, height)
            if unkc:
                if nc == 1:
                    mode = 'L'
                elif nc == 2:
                    mode = 'LA'
                elif nc == 3:
                    mode = 'RGB'
                elif nc == 4:
                    mode = 'RGBA'
                break
        elif tbox == b'colr':
            meth, prec, approx = struct.unpack('>BBB', content[:3])
            if meth == 1:
                cs = struct.unpack('>I', content[3:7])[0]
                if cs == 16:   # sRGB
                    if nc == 3:
                        mode = 'RGB'
                    elif nc == 4:
                        mode = 'RGBA'
                    break
                elif cs == 17: # grayscale
                    if nc == 1:
                        mode = 'L'
                    elif nc == 2:
                        mode = 'LA'
                    break
                elif cs == 18: # sYCC
                    if nc == 3:
                        mode = 'RGB'
                    elif nc == 4:
                        mode == 'RGBA'
                    break

    return (size, mode)

##
# Image plugin for JPEG2000 images.

class Jpeg2KImageFile(ImageFile.ImageFile):
    format = "JPEG2000"
    format_description = "JPEG 2000 (ISO 15444)"

    def _open(self):
        sig = self.fp.read(4)
        if sig == b'\xff\x4f\xff\x51':
            self.codec = "j2k"
            try:
                self.size, self.mode = _parse_codestream(self.fp)
            except Exception as e:
                print e
        else:
            sig = sig + self.fp.read(8)
        
            if sig == b'\x00\x00\x00\x0cjP  \x0d\x0a\x87\x0a':
                self.codec = "jp2"
                self.size, self.mode = _parse_jp2_header(self.fp)
            else:
                raise SyntaxError('not a JPEG 2000 file')
        
        if self.size is None or self.mode is None:
            raise SyntaxError('unable to determine size/mode')
        
        self.reduce = 0
        self.layers = 0

        self.tile = [('jpeg2k', (0, 0) + self.size, 0,
                      (self.codec, self.reduce, self.layers))]

    def load(self):
        if self.reduce:
            power = 1 << self.reduce
            adjust = power >> 1
            self.size = ((self.size[0] + adjust) / power,
                         (self.size[1] + adjust) / power)
        
        ImageFile.ImageFile.load(self)
        
def _accept(prefix):
    return (prefix[:4] == b'\xff\x4f\xff\x51'
            or prefix[:12] == b'\x00\x00\x00\x0cjP  \x0d\x0a\x87\x0a')

# ------------------------------------------------------------
# Save support

def _save(im, fp, filename):
    if filename.endswith('.j2k'):
        kind = 'j2k'
    else:
        kind = 'jp2'
        
    ImageFile._save(im, fp, [('jpeg2k', (0, 0)+im.size, 0, kind)])
    
# ------------------------------------------------------------
# Registry stuff

Image.register_open('JPEG2000', Jpeg2KImageFile, _accept)
Image.register_save('JPEG2000', _save)

Image.register_extension('JPEG2000', '.jp2')
Image.register_extension('JPEG2000', '.j2k')
Image.register_extension('JPEG2000', '.jpc')
Image.register_extension('JPEG2000', '.jpf')
Image.register_extension('JPEG2000', '.jpx')
Image.register_extension('JPEG2000', '.j2c')

Image.register_mime('JPEG2000', 'image/jp2')
Image.register_mime('JPEG2000', 'image/jpx')
