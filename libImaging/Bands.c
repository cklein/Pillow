/*
 * The Python Imaging Library
 * $Id$
 *
 * stuff to extract and paste back individual bands
 *
 * history:
 * 1996-03-20 fl   Created
 * 1997-08-27 fl   Fixed putband for single band targets.
 * 2003-09-26 fl   Fixed getband/putband for 2-band images (LA, PA).
 *
 * Copyright (c) 1997-2003 by Secret Labs AB.
 * Copyright (c) 1996-1997 by Fredrik Lundh.
 *
 * See the README file for details on usage and redistribution.
 */


#include "Imaging.h"


#define CLIP(x) ((x) <= 0 ? 0 : (x) < 256 ? (x) : 255)


#ifdef WORDS_BIGENDIAN
    #define MAKE_UINT32(u0, u1, u2, u3) (u3 | (u2<<8) | (u1<<16) | (u0<<24))
#else
    #define MAKE_UINT32(u0, u1, u2, u3) (u0 | (u1<<8) | (u2<<16) | (u3<<24))
#endif


Imaging
ImagingGetBand(Imaging imIn, int band)
{
    Imaging imOut;
    int x, y;

    /* Check arguments */
    if (!imIn || imIn->type != IMAGING_TYPE_UINT8)
        return (Imaging) ImagingError_ModeError();

    if (band < 0 || band >= imIn->bands)
        return (Imaging) ImagingError_ValueError("band index out of range");

    /* Shortcuts */
    if (imIn->bands == 1)
        return ImagingCopy(imIn);

    /* Special case for LXXA etc */
    if (imIn->bands == 2 && band == 1)
        band = 3;

    imOut = ImagingNew("L", imIn->xsize, imIn->ysize);
    if (!imOut)
        return NULL;

    /* Extract band from image */
    for (y = 0; y < imIn->ysize; y++) {
        UINT8* in = (UINT8*) imIn->image[y] + band;
        UINT8* out = imOut->image8[y];
        for (x = 0; x < imIn->xsize; x++) {
            out[x] = *in;
            in += 4;
        }
    }

    return imOut;
}

Imaging
ImagingPutBand(Imaging imOut, Imaging imIn, int band)
{
    int x, y;

    /* Check arguments */
    if (!imIn || imIn->bands != 1 || !imOut)
        return (Imaging) ImagingError_ModeError();

    if (band < 0 || band >= imOut->bands)
        return (Imaging) ImagingError_ValueError("band index out of range");

    if (imIn->type  != imOut->type  ||
        imIn->xsize != imOut->xsize ||
        imIn->ysize != imOut->ysize)
        return (Imaging) ImagingError_Mismatch();

    /* Shortcuts */
    if (imOut->bands == 1)
        return ImagingCopy2(imOut, imIn);

    /* Special case for LXXA etc */
    if (imOut->bands == 2 && band == 1)
        band = 3;

    /* Insert band into image */
    for (y = 0; y < imIn->ysize; y++) {
        UINT8* in = imIn->image8[y];
        UINT8* out = (UINT8*) imOut->image[y] + band;
        for (x = 0; x < imIn->xsize; x++) {
            *out = in[x];
            out += 4;
        }
    }

    return imOut;
}

Imaging
ImagingFillBand(Imaging imOut, int band, int color)
{
    int x, y;

    /* Check arguments */
    if (!imOut || imOut->type != IMAGING_TYPE_UINT8)
        return (Imaging) ImagingError_ModeError();

    if (band < 0 || band >= imOut->bands)
        return (Imaging) ImagingError_ValueError("band index out of range");

    /* Special case for LXXA etc */
    if (imOut->bands == 2 && band == 1)
        band = 3;

    color = CLIP(color);

    /* Insert color into image */
    for (y = 0; y < imOut->ysize; y++) {
        UINT8* out = (UINT8*) imOut->image[y] + band;
        for (x = 0; x < imOut->xsize; x++) {
            *out = (UINT8) color;
            out += 4;
        }
    }

    return imOut;
}

Imaging
ImagingMerge(const char* mode, Imaging bands[4])
{
    int i, x, y;
    int bandsCount = 0;
    Imaging imOut;
    Imaging firstBand;

    firstBand = bands[0];
    if ( ! firstBand) {
        return (Imaging) ImagingError_ValueError("wrong number of bands");
    }

    for (i = 0; i < 4; ++i) {
        if ( ! bands[i]) {
            break;
        }
        if (bands[i]->bands != 1) {
            return (Imaging) ImagingError_ModeError();
        }
        if (bands[i]->xsize != firstBand->xsize
            || bands[i]->ysize != firstBand->ysize) {
            return (Imaging) ImagingError_Mismatch();
        }
    }
    bandsCount = i;

    imOut = ImagingNew(mode, firstBand->xsize, firstBand->ysize);
    if ( ! imOut)
        return NULL;

    if (imOut->bands != bandsCount) {
        ImagingDelete(imOut);
        return (Imaging) ImagingError_ValueError("wrong number of bands");
    }

    if (imOut->bands == 1)
        return ImagingCopy2(imOut, firstBand);

    if (imOut->bands == 2) {
        for (y = 0; y < imOut->ysize; y++) {
            UINT8* in0 = bands[0]->image8[y];
            UINT8* in1 = bands[1]->image8[y];
            UINT32* out = (UINT32*) imOut->image32[y];
            for (x = 0; x < imOut->xsize; x++) {
                out[x] = MAKE_UINT32(in0[x], 0, 0, in1[x]);
            }
        }
    } else if (imOut->bands == 3) {
        for (y = 0; y < imOut->ysize; y++) {
            UINT8* in0 = bands[0]->image8[y];
            UINT8* in1 = bands[1]->image8[y];
            UINT8* in2 = bands[2]->image8[y];
            UINT32* out = (UINT32*) imOut->image32[y];
            for (x = 0; x < imOut->xsize; x++) {
                out[x] = MAKE_UINT32(in0[x], in1[x], in2[x], 0);
            }
        }
    } else if (imOut->bands == 4) {
        for (y = 0; y < imOut->ysize; y++) {
            UINT8* in0 = bands[0]->image8[y];
            UINT8* in1 = bands[1]->image8[y];
            UINT8* in2 = bands[2]->image8[y];
            UINT8* in3 = bands[3]->image8[y];
            UINT32* out = (UINT32*) imOut->image32[y];
            for (x = 0; x < imOut->xsize; x++) {
                out[x] = MAKE_UINT32(in0[x], in1[x], in2[x], in3[x]);
            }
        }
    }

    return imOut;
}
