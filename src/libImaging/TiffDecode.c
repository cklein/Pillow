/*
 * The Python Imaging Library.
 * $Id: //modules/pil/libImaging/TiffDecode.c#1 $
 *
 * LibTiff-based Group3 and Group4 decoder
 *
 *
 * started modding to use non-private tiff functions to port to libtiff 4.x
 * eds 3/12/12
 *
 */

#include "Imaging.h"

#ifdef HAVE_LIBTIFF

#ifndef uint
#define uint uint32
#endif

#include "TiffDecode.h"

/* Convert C file descriptor to WinApi HFILE if LibTiff was compiled with tif_win32.c
 *
 * This cast is safe, as the top 32-bits of HFILE are guaranteed to be zero,
 * see
 * https://docs.microsoft.com/en-us/windows/win32/winprog64/interprocess-communication
 */
#ifndef USE_WIN32_FILEIO
#define fd_to_tiff_fd(fd) (fd)
#else
#define fd_to_tiff_fd(fd) ((int)_get_osfhandle(fd))
#endif

void
dump_state(const TIFFSTATE *state) {
    TRACE(
        ("State: Location %u size %d eof %d data: %p ifd: %d\n",
         (uint)state->loc,
         (int)state->size,
         (uint)state->eof,
         state->data,
         state->ifd));
}

/*
  procs for TIFFOpenClient
*/

tsize_t
_tiffReadProc(thandle_t hdata, tdata_t buf, tsize_t size) {
    TIFFSTATE *state = (TIFFSTATE *)hdata;
    tsize_t to_read;

    TRACE(("_tiffReadProc: %d \n", (int)size));
    dump_state(state);

    if (state->loc > state->eof) {
        TIFFError("_tiffReadProc", "Invalid Read at loc %d, eof: %d", state->loc, state->eof);
        return 0;
    }
    to_read = min(size, min(state->size, (tsize_t)state->eof) - (tsize_t)state->loc);
    TRACE(("to_read: %d\n", (int)to_read));

    _TIFFmemcpy(buf, (UINT8 *)state->data + state->loc, to_read);
    state->loc += (toff_t)to_read;

    TRACE(("location: %u\n", (uint)state->loc));
    return to_read;
}

tsize_t
_tiffWriteProc(thandle_t hdata, tdata_t buf, tsize_t size) {
    TIFFSTATE *state = (TIFFSTATE *)hdata;
    tsize_t to_write;

    TRACE(("_tiffWriteProc: %d \n", (int)size));
    dump_state(state);

    to_write = min(size, state->size - (tsize_t)state->loc);
    if (state->flrealloc && size > to_write) {
        tdata_t new_data;
        tsize_t newsize = state->size;
        while (newsize < (size + state->size)) {
            if (newsize > INT_MAX - 64 * 1024) {
                return 0;
            }
            newsize += 64 * 1024;
            // newsize*=2; // UNDONE, by 64k chunks?
        }
        TRACE(("Reallocing in write to %d bytes\n", (int)newsize));
        /* malloc check ok, overflow checked above */
        new_data = realloc(state->data, newsize);
        if (!new_data) {
            // fail out
            return 0;
        }
        state->data = new_data;
        state->size = newsize;
        to_write = size;
    }

    TRACE(("to_write: %d\n", (int)to_write));

    _TIFFmemcpy((UINT8 *)state->data + state->loc, buf, to_write);
    state->loc += (toff_t)to_write;
    state->eof = max(state->loc, state->eof);

    dump_state(state);
    return to_write;
}

toff_t
_tiffSeekProc(thandle_t hdata, toff_t off, int whence) {
    TIFFSTATE *state = (TIFFSTATE *)hdata;

    TRACE(("_tiffSeekProc: off: %u whence: %d \n", (uint)off, whence));
    dump_state(state);
    switch (whence) {
        case 0:
            state->loc = off;
            break;
        case 1:
            state->loc += off;
            break;
        case 2:
            state->loc = state->eof + off;
            break;
    }
    dump_state(state);
    return state->loc;
}

int
_tiffCloseProc(thandle_t hdata) {
    TIFFSTATE *state = (TIFFSTATE *)hdata;

    TRACE(("_tiffCloseProc \n"));
    dump_state(state);

    return 0;
}

toff_t
_tiffSizeProc(thandle_t hdata) {
    TIFFSTATE *state = (TIFFSTATE *)hdata;

    TRACE(("_tiffSizeProc \n"));
    dump_state(state);

    return (toff_t)state->size;
}

int
_tiffMapProc(thandle_t hdata, tdata_t *pbase, toff_t *psize) {
    TIFFSTATE *state = (TIFFSTATE *)hdata;

    TRACE(("_tiffMapProc input size: %u, data: %p\n", (uint)*psize, *pbase));
    dump_state(state);

    *pbase = state->data;
    *psize = state->size;
    TRACE(("_tiffMapProc returning size: %u, data: %p\n", (uint)*psize, *pbase));
    return (1);
}

int
_tiffNullMapProc(thandle_t hdata, tdata_t *pbase, toff_t *psize) {
    (void)hdata;
    (void)pbase;
    (void)psize;
    return (0);
}

void
_tiffUnmapProc(thandle_t hdata, tdata_t base, toff_t size) {
    TRACE(("_tiffUnMapProc\n"));
    (void)hdata;
    (void)base;
    (void)size;
}

int
ImagingLibTiffInit(ImagingCodecState state, int fp, uint32 offset) {
    TIFFSTATE *clientstate = (TIFFSTATE *)state->context;

    TRACE(("initing libtiff\n"));
    TRACE(("filepointer: %d \n", fp));
    TRACE(
        ("State: count %d, state %d, x %d, y %d, ystep %d\n",
         state->count,
         state->state,
         state->x,
         state->y,
         state->ystep));
    TRACE(
        ("State: xsize %d, ysize %d, xoff %d, yoff %d \n",
         state->xsize,
         state->ysize,
         state->xoff,
         state->yoff));
    TRACE(("State: bits %d, bytes %d \n", state->bits, state->bytes));
    TRACE(("State: context %p \n", state->context));

    clientstate->loc = 0;
    clientstate->size = 0;
    clientstate->data = 0;
    clientstate->fp = fp;
    clientstate->ifd = offset;
    clientstate->eof = 0;

    return 1;
}

int
_decodeStripYCbCr(Imaging im, ImagingCodecState state, TIFF *tiff) {
    // To avoid dealing with YCbCr subsampling, let libtiff handle it
    // Use a TIFFRGBAImage wrapping the tiff image, and let libtiff handle
    // all of the conversion. Metadata read from the TIFFRGBAImage could
    // be different from the metadata that the base tiff returns.

    INT32 strip_row;
    UINT8 *new_data;
    UINT32 rows_per_strip, row_byte_size, rows_to_read;
    int ret;
    TIFFRGBAImage img;
    char emsg[1024] = "";

    ret = TIFFGetFieldDefaulted(tiff, TIFFTAG_ROWSPERSTRIP, &rows_per_strip);
    if (ret != 1) {
        rows_per_strip = state->ysize;
    }
    TRACE(("RowsPerStrip: %u \n", rows_per_strip));

    if (!(TIFFRGBAImageOK(tiff, emsg) && TIFFRGBAImageBegin(&img, tiff, 0, emsg))) {
        TRACE(("Decode error, msg: %s", emsg));
        state->errcode = IMAGING_CODEC_BROKEN;
        // nothing to clean up, just return
        return -1;
    }

    img.req_orientation = ORIENTATION_TOPLEFT;
    img.col_offset = 0;

    if (state->xsize != img.width || state->ysize != img.height) {
        TRACE(
            ("Inconsistent Image Error: %d =? %d, %d =? %d",
             state->xsize,
             img.width,
             state->ysize,
             img.height));
        state->errcode = IMAGING_CODEC_BROKEN;
        goto decodeycbcr_err;
    }

    /* overflow check for row byte size */
    if (INT_MAX / 4 < img.width) {
        state->errcode = IMAGING_CODEC_MEMORY;
        goto decodeycbcr_err;
    }

    // TiffRGBAImages are 32bits/pixel.
    row_byte_size = img.width * 4;

    /* overflow check for realloc */
    if (INT_MAX / row_byte_size < rows_per_strip) {
        state->errcode = IMAGING_CODEC_MEMORY;
        goto decodeycbcr_err;
    }

    state->bytes = rows_per_strip * row_byte_size;

    TRACE(("StripSize: %d \n", state->bytes));

    /* realloc to fit whole strip */
    /* malloc check above */
    new_data = realloc(state->buffer, state->bytes);
    if (!new_data) {
        state->errcode = IMAGING_CODEC_MEMORY;
        goto decodeycbcr_err;
    }

    state->buffer = new_data;

    for (; state->y < state->ysize; state->y += rows_per_strip) {
        img.row_offset = state->y;
        rows_to_read = min(rows_per_strip, img.height - state->y);

        if (!TIFFRGBAImageGet(&img, (UINT32 *)state->buffer, img.width, rows_to_read)) {
            TRACE(("Decode Error, y: %d\n", state->y));
            state->errcode = IMAGING_CODEC_BROKEN;
            goto decodeycbcr_err;
        }

        TRACE(("Decoded strip for row %d \n", state->y));

        // iterate over each row in the strip and stuff data into image
        for (strip_row = 0;
             strip_row < min((INT32)rows_per_strip, state->ysize - state->y);
             strip_row++) {
            TRACE(("Writing data into line %d ; \n", state->y + strip_row));

            // UINT8 * bbb = state->buffer + strip_row * (state->bytes /
            // rows_per_strip); TRACE(("chars: %x %x %x %x\n", ((UINT8 *)bbb)[0],
            // ((UINT8 *)bbb)[1], ((UINT8 *)bbb)[2], ((UINT8 *)bbb)[3]));

            state->shuffle(
                (UINT8 *)im->image[state->y + state->yoff + strip_row] +
                    state->xoff * im->pixelsize,
                state->buffer + strip_row * row_byte_size,
                state->xsize);
        }
    }

decodeycbcr_err:
    TIFFRGBAImageEnd(&img);
    if (state->errcode != 0) {
        return -1;
    }
    return 0;
}

int
_decodeStrip(Imaging im, ImagingCodecState state, TIFF *tiff) {
    INT32 strip_row;
    UINT8 *new_data;
    UINT32 rows_per_strip, row_byte_size;
    int ret;

    ret = TIFFGetField(tiff, TIFFTAG_ROWSPERSTRIP, &rows_per_strip);
    if (ret != 1) {
        rows_per_strip = state->ysize;
    }
    TRACE(("RowsPerStrip: %u \n", rows_per_strip));

    // We could use TIFFStripSize, but for YCbCr data it returns subsampled data size
    row_byte_size = (state->xsize * state->bits + 7) / 8;

    /* overflow check for realloc */
    if (INT_MAX / row_byte_size < rows_per_strip) {
        state->errcode = IMAGING_CODEC_MEMORY;
        return -1;
    }

    state->bytes = rows_per_strip * row_byte_size;

    TRACE(("StripSize: %d \n", state->bytes));

    if (TIFFStripSize(tiff) > state->bytes) {
        // If the strip size as expected by LibTiff isn't what we're expecting, abort.
        // man:   TIFFStripSize returns the equivalent size for a strip of data as it
        // would be returned in a
        //        call to TIFFReadEncodedStrip ...

        state->errcode = IMAGING_CODEC_MEMORY;
        return -1;
    }

    /* realloc to fit whole strip */
    /* malloc check above */
    new_data = realloc(state->buffer, state->bytes);
    if (!new_data) {
        state->errcode = IMAGING_CODEC_MEMORY;
        return -1;
    }

    state->buffer = new_data;

    for (; state->y < state->ysize; state->y += rows_per_strip) {
        if (TIFFReadEncodedStrip(
                tiff,
                TIFFComputeStrip(tiff, state->y, 0),
                (tdata_t)state->buffer,
                -1) == -1) {
            TRACE(("Decode Error, strip %d\n", TIFFComputeStrip(tiff, state->y, 0)));
            state->errcode = IMAGING_CODEC_BROKEN;
            return -1;
        }

        TRACE(("Decoded strip for row %d \n", state->y));

        // iterate over each row in the strip and stuff data into image
        for (strip_row = 0;
             strip_row < min((INT32)rows_per_strip, state->ysize - state->y);
             strip_row++) {
            TRACE(("Writing data into line %d ; \n", state->y + strip_row));

            // UINT8 * bbb = state->buffer + strip_row * (state->bytes /
            // rows_per_strip); TRACE(("chars: %x %x %x %x\n", ((UINT8 *)bbb)[0],
            // ((UINT8 *)bbb)[1], ((UINT8 *)bbb)[2], ((UINT8 *)bbb)[3]));

            state->shuffle(
                (UINT8 *)im->image[state->y + state->yoff + strip_row] +
                    state->xoff * im->pixelsize,
                state->buffer + strip_row * row_byte_size,
                state->xsize);
        }
    }
    return 0;
}

int
ImagingLibTiffDecode(
    Imaging im, ImagingCodecState state, UINT8 *buffer, Py_ssize_t bytes) {
    TIFFSTATE *clientstate = (TIFFSTATE *)state->context;
    char *filename = "tempfile.tif";
    char *mode = "r";
    TIFF *tiff;
    uint16 photometric = 0;  // init to not PHOTOMETRIC_YCBCR
    int isYCbCr = 0;

    /* buffer is the encoded file, bytes is the length of the encoded file */
    /*     it all ends up in state->buffer, which is a uint8* from Imaging.h */

    TRACE(("in decoder: bytes %d\n", bytes));
    TRACE(
        ("State: count %d, state %d, x %d, y %d, ystep %d\n",
         state->count,
         state->state,
         state->x,
         state->y,
         state->ystep));
    TRACE(
        ("State: xsize %d, ysize %d, xoff %d, yoff %d \n",
         state->xsize,
         state->ysize,
         state->xoff,
         state->yoff));
    TRACE(("State: bits %d, bytes %d \n", state->bits, state->bytes));
    TRACE(
        ("Buffer: %p: %c%c%c%c\n",
         buffer,
         (char)buffer[0],
         (char)buffer[1],
         (char)buffer[2],
         (char)buffer[3]));
    TRACE(
        ("State->Buffer: %c%c%c%c\n",
         (char)state->buffer[0],
         (char)state->buffer[1],
         (char)state->buffer[2],
         (char)state->buffer[3]));
    TRACE(
        ("Image: mode %s, type %d, bands: %d, xsize %d, ysize %d \n",
         im->mode,
         im->type,
         im->bands,
         im->xsize,
         im->ysize));
    TRACE(
        ("Image: image8 %p, image32 %p, image %p, block %p \n",
         im->image8,
         im->image32,
         im->image,
         im->block));
    TRACE(("Image: pixelsize: %d, linesize %d \n", im->pixelsize, im->linesize));

    dump_state(clientstate);
    clientstate->size = bytes;
    clientstate->eof = clientstate->size;
    clientstate->loc = 0;
    clientstate->data = (tdata_t)buffer;
    clientstate->flrealloc = 0;
    dump_state(clientstate);

    TIFFSetWarningHandler(NULL);
    TIFFSetWarningHandlerExt(NULL);

    if (clientstate->fp) {
        TRACE(("Opening using fd: %d\n", clientstate->fp));
        lseek(clientstate->fp, 0, SEEK_SET);  // Sometimes, I get it set to the end.
        tiff = TIFFFdOpen(fd_to_tiff_fd(clientstate->fp), filename, mode);
    } else {
        TRACE(("Opening from string\n"));
        tiff = TIFFClientOpen(
            filename,
            mode,
            (thandle_t)clientstate,
            _tiffReadProc,
            _tiffWriteProc,
            _tiffSeekProc,
            _tiffCloseProc,
            _tiffSizeProc,
            _tiffMapProc,
            _tiffUnmapProc);
    }

    if (!tiff) {
        TRACE(("Error, didn't get the tiff\n"));
        state->errcode = IMAGING_CODEC_BROKEN;
        return -1;
    }

    if (clientstate->ifd) {
        int rv;
        uint32 ifdoffset = clientstate->ifd;
        TRACE(("reading tiff ifd %u\n", ifdoffset));
        rv = TIFFSetSubDirectory(tiff, ifdoffset);
        if (!rv) {
            TRACE(("error in TIFFSetSubDirectory"));
            goto decode_err;
        }
    }

    TIFFGetField(tiff, TIFFTAG_PHOTOMETRIC, &photometric);
    isYCbCr = photometric == PHOTOMETRIC_YCBCR;

    if (TIFFIsTiled(tiff)) {
        INT32 x, y, tile_y;
        UINT32 tile_width, tile_length, current_tile_length, current_line,
            current_tile_width, row_byte_size;
        UINT8 *new_data;

        TIFFGetField(tiff, TIFFTAG_TILEWIDTH, &tile_width);
        TIFFGetField(tiff, TIFFTAG_TILELENGTH, &tile_length);

        /* overflow check for row_byte_size calculation */
        if ((UINT32)INT_MAX / state->bits < tile_width) {
            state->errcode = IMAGING_CODEC_MEMORY;
            goto decode_err;
        }

        if (isYCbCr) {
            row_byte_size = tile_width * 4;
            /* sanity check, we use this value in shuffle below */
            if (im->pixelsize != 4) {
                state->errcode = IMAGING_CODEC_BROKEN;
                goto decode_err;
            }
        } else {
            // We could use TIFFTileSize, but for YCbCr data it returns subsampled data
            // size
            row_byte_size = (tile_width * state->bits + 7) / 8;
        }

        /* overflow check for realloc */
        if (INT_MAX / row_byte_size < tile_length) {
            state->errcode = IMAGING_CODEC_MEMORY;
            goto decode_err;
        }

        state->bytes = row_byte_size * tile_length;

        if (TIFFTileSize(tiff) > state->bytes) {
            // If the strip size as expected by LibTiff isn't what we're expecting,
            // abort.
            state->errcode = IMAGING_CODEC_MEMORY;
            goto decode_err;
        }

        /* realloc to fit whole tile */
        /* malloc check above */
        new_data = realloc(state->buffer, state->bytes);
        if (!new_data) {
            state->errcode = IMAGING_CODEC_MEMORY;
            goto decode_err;
        }

        state->buffer = new_data;

        TRACE(("TIFFTileSize: %d\n", state->bytes));

        for (y = state->yoff; y < state->ysize; y += tile_length) {
            for (x = state->xoff; x < state->xsize; x += tile_width) {
                /* Sanity Check. Apparently in some cases, the TiffReadRGBA* functions
                   have a different view of the size of the tiff than we're getting from
                   other functions. So, we need to check here. 
                */
                if (!TIFFCheckTile(tiff, x, y, 0, 0)) {
                    TRACE(("Check Tile Error, Tile at %dx%d\n", x, y));
                    state->errcode = IMAGING_CODEC_BROKEN;
                    goto decode_err;
                }
                if (isYCbCr) {
                    /* To avoid dealing with YCbCr subsampling, let libtiff handle it */
                    if (!TIFFReadRGBATile(tiff, x, y, (UINT32 *)state->buffer)) {
                        TRACE(("Decode Error, Tile at %dx%d\n", x, y));
                        state->errcode = IMAGING_CODEC_BROKEN;
                        goto decode_err;
                    }
                } else {
                    if (TIFFReadTile(tiff, (tdata_t)state->buffer, x, y, 0, 0) == -1) {
                        TRACE(("Decode Error, Tile at %dx%d\n", x, y));
                        state->errcode = IMAGING_CODEC_BROKEN;
                        goto decode_err;
                    }
                }

                TRACE(("Read tile at %dx%d; \n\n", x, y));

                current_tile_width = min((INT32)tile_width, state->xsize - x);
                current_tile_length = min((INT32)tile_length, state->ysize - y);
                // iterate over each line in the tile and stuff data into image
                for (tile_y = 0; tile_y < current_tile_length; tile_y++) {
                    TRACE(
                        ("Writing tile data at %dx%d using tile_width: %d; \n",
                         tile_y + y,
                         x,
                         current_tile_width));

                    // UINT8 * bbb = state->buffer + tile_y * row_byte_size;
                    // TRACE(("chars: %x%x%x%x\n", ((UINT8 *)bbb)[0], ((UINT8 *)bbb)[1],
                    // ((UINT8 *)bbb)[2], ((UINT8 *)bbb)[3]));
                    /*
                     * For some reason the TIFFReadRGBATile() function
                     * chooses the lower left corner as the origin.
                     * Vertically mirror by shuffling the scanlines
                     * backwards
                     */

                    if (isYCbCr) {
                        current_line = tile_length - tile_y - 1;
                    } else {
                        current_line = tile_y;
                    }

                    state->shuffle(
                        (UINT8 *)im->image[tile_y + y] + x * im->pixelsize,
                        state->buffer + current_line * row_byte_size,
                        current_tile_width);
                }
            }
        }
    } else {
        if (!isYCbCr) {
            _decodeStrip(im, state, tiff);
        } else {
            _decodeStripYCbCr(im, state, tiff);
        }
    }

decode_err:
    TIFFClose(tiff);
    TRACE(("Done Decoding, Returning \n"));
    // Returning -1 here to force ImageFile.load to break, rather than
    // even think about looping back around.
    return -1;
}

int
ImagingLibTiffEncodeInit(ImagingCodecState state, char *filename, int fp) {
    // Open the FD or the pointer as a tiff file, for writing.
    // We may have to do some monkeying around to make this really work.
    // If we have a fp, then we're good.
    // If we have a memory string, we're probably going to have to malloc, then
    // shuffle bytes into the writescanline process.
    // Going to have to deal with the directory as well.

    TIFFSTATE *clientstate = (TIFFSTATE *)state->context;
    int bufsize = 64 * 1024;
    char *mode = "w";

    TRACE(("initing libtiff\n"));
    TRACE(("Filename %s, filepointer: %d \n", filename, fp));
    TRACE(
        ("State: count %d, state %d, x %d, y %d, ystep %d\n",
         state->count,
         state->state,
         state->x,
         state->y,
         state->ystep));
    TRACE(
        ("State: xsize %d, ysize %d, xoff %d, yoff %d \n",
         state->xsize,
         state->ysize,
         state->xoff,
         state->yoff));
    TRACE(("State: bits %d, bytes %d \n", state->bits, state->bytes));
    TRACE(("State: context %p \n", state->context));

    clientstate->loc = 0;
    clientstate->size = 0;
    clientstate->eof = 0;
    clientstate->data = 0;
    clientstate->flrealloc = 0;
    clientstate->fp = fp;

    state->state = 0;

    if (fp) {
        TRACE(("Opening using fd: %d for writing \n", clientstate->fp));
        clientstate->tiff = TIFFFdOpen(fd_to_tiff_fd(clientstate->fp), filename, mode);
    } else {
        // malloc a buffer to write the tif, we're going to need to realloc or something
        // if we need bigger.
        TRACE(("Opening a buffer for writing \n"));
        /* malloc check ok, small constant allocation */
        clientstate->data = malloc(bufsize);
        clientstate->size = bufsize;
        clientstate->flrealloc = 1;

        if (!clientstate->data) {
            TRACE(("Error, couldn't allocate a buffer of size %d\n", bufsize));
            return 0;
        }

        clientstate->tiff = TIFFClientOpen(
            filename,
            mode,
            (thandle_t)clientstate,
            _tiffReadProc,
            _tiffWriteProc,
            _tiffSeekProc,
            _tiffCloseProc,
            _tiffSizeProc,
            _tiffNullMapProc,
            _tiffUnmapProc); /*force no mmap*/
    }

    if (!clientstate->tiff) {
        TRACE(("Error, couldn't open tiff file\n"));
        return 0;
    }

    return 1;
}

int
ImagingLibTiffMergeFieldInfo(
    ImagingCodecState state, TIFFDataType field_type, int key, int is_var_length) {
    // Refer to libtiff docs (http://www.simplesystems.org/libtiff/addingtags.html)
    TIFFSTATE *clientstate = (TIFFSTATE *)state->context;
    uint32 n;
    int status = 0;

    // custom fields added with ImagingLibTiffMergeFieldInfo are only used for
    // decoding, ignore readcount;
    int readcount = 1;
    // we support writing a single value, or a variable number of values
    int writecount = 1;
    // whether the first value should encode the number of values.
    int passcount = 0;

    TIFFFieldInfo info[] = {
        {key,
         readcount,
         writecount,
         field_type,
         FIELD_CUSTOM,
         1,
         passcount,
         "CustomField"}};

    if (is_var_length) {
        info[0].field_writecount = -1;
    }

    if (is_var_length && field_type != TIFF_ASCII) {
        info[0].field_passcount = 1;
    }

    n = sizeof(info) / sizeof(info[0]);

    // Test for libtiff 4.0 or later, excluding libtiff 3.9.6 and 3.9.7
#if TIFFLIB_VERSION >= 20111221 && TIFFLIB_VERSION != 20120218 && \
    TIFFLIB_VERSION != 20120922
    status = TIFFMergeFieldInfo(clientstate->tiff, info, n);
#else
    TIFFMergeFieldInfo(clientstate->tiff, info, n);
#endif
    return status;
}

int
ImagingLibTiffSetField(ImagingCodecState state, ttag_t tag, ...) {
    // after tif_dir.c->TIFFSetField.
    TIFFSTATE *clientstate = (TIFFSTATE *)state->context;
    va_list ap;
    int status;

    va_start(ap, tag);
    status = TIFFVSetField(clientstate->tiff, tag, ap);
    va_end(ap);
    return status;
}

int
ImagingLibTiffEncode(Imaging im, ImagingCodecState state, UINT8 *buffer, int bytes) {
    /* One shot encoder. Encode everything to the tiff in the clientstate.
       If we're running off of a FD, then run once, we're good, everything
       ends up in the file, we close and we're done.

       If we're going to memory, then we need to write the whole file into memory, then
       parcel it back out to the pystring buffer bytes at a time.

    */

    TIFFSTATE *clientstate = (TIFFSTATE *)state->context;
    TIFF *tiff = clientstate->tiff;

    TRACE(("in encoder: bytes %d\n", bytes));
    TRACE(
        ("State: count %d, state %d, x %d, y %d, ystep %d\n",
         state->count,
         state->state,
         state->x,
         state->y,
         state->ystep));
    TRACE(
        ("State: xsize %d, ysize %d, xoff %d, yoff %d \n",
         state->xsize,
         state->ysize,
         state->xoff,
         state->yoff));
    TRACE(("State: bits %d, bytes %d \n", state->bits, state->bytes));
    TRACE(
        ("Buffer: %p: %c%c%c%c\n",
         buffer,
         (char)buffer[0],
         (char)buffer[1],
         (char)buffer[2],
         (char)buffer[3]));
    TRACE(
        ("State->Buffer: %c%c%c%c\n",
         (char)state->buffer[0],
         (char)state->buffer[1],
         (char)state->buffer[2],
         (char)state->buffer[3]));
    TRACE(
        ("Image: mode %s, type %d, bands: %d, xsize %d, ysize %d \n",
         im->mode,
         im->type,
         im->bands,
         im->xsize,
         im->ysize));
    TRACE(
        ("Image: image8 %p, image32 %p, image %p, block %p \n",
         im->image8,
         im->image32,
         im->image,
         im->block));
    TRACE(("Image: pixelsize: %d, linesize %d \n", im->pixelsize, im->linesize));

    dump_state(clientstate);

    if (state->state == 0) {
        TRACE(("Encoding line bt line"));
        while (state->y < state->ysize) {
            state->shuffle(
                state->buffer,
                (UINT8 *)im->image[state->y + state->yoff] +
                    state->xoff * im->pixelsize,
                state->xsize);

            if (TIFFWriteScanline(
                    tiff, (tdata_t)(state->buffer), (uint32)state->y, 0) == -1) {
                TRACE(("Encode Error, row %d\n", state->y));
                state->errcode = IMAGING_CODEC_BROKEN;
                TIFFClose(tiff);
                if (!clientstate->fp) {
                    free(clientstate->data);
                }
                return -1;
            }
            state->y++;
        }

        if (state->y == state->ysize) {
            state->state = 1;

            TRACE(("Flushing \n"));
            if (!TIFFFlush(tiff)) {
                TRACE(("Error flushing the tiff"));
                // likely reason is memory.
                state->errcode = IMAGING_CODEC_MEMORY;
                TIFFClose(tiff);
                if (!clientstate->fp) {
                    free(clientstate->data);
                }
                return -1;
            }
            TRACE(("Closing \n"));
            TIFFClose(tiff);
            // reset the clientstate metadata to use it to read out the buffer.
            clientstate->loc = 0;
            clientstate->size = clientstate->eof;  // redundant?
        }
    }

    if (state->state == 1 && !clientstate->fp) {
        int read = (int)_tiffReadProc(clientstate, (tdata_t)buffer, (tsize_t)bytes);
        TRACE(
            ("Buffer: %p: %c%c%c%c\n",
             buffer,
             (char)buffer[0],
             (char)buffer[1],
             (char)buffer[2],
             (char)buffer[3]));
        if (clientstate->loc == clientstate->eof) {
            TRACE(("Hit EOF, calling an end, freeing data"));
            state->errcode = IMAGING_CODEC_END;
            free(clientstate->data);
        }
        return read;
    }

    state->errcode = IMAGING_CODEC_END;
    return 0;
}

const char *
ImagingTiffVersion(void) {
    return TIFFGetVersion();
}

#endif
