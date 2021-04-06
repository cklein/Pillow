/*
 * The Python Imaging Library
 *
 * a simple math add-on for the Python Imaging Library
 *
 * history:
 * 1999-02-15 fl   Created
 * 2005-05-05 fl   Simplified and cleaned up for PIL 1.1.6
 *
 * Copyright (c) 1999-2005 by Secret Labs AB
 * Copyright (c) 2005 by Fredrik Lundh
 *
 * See the README file for information on usage and redistribution.
 */

#include "hpy.h"

#include "libImaging/Imaging.h"

#include "math.h"
#include "float.h"

#define MAX_INT32 2147483647.0
#define MIN_INT32 -2147483648.0

#define UNOP(name, op, type)                   \
    void name(Imaging out, Imaging im1) {      \
        int x, y;                              \
        for (y = 0; y < out->ysize; y++) {     \
            type *p0 = (type *)out->image[y];  \
            type *p1 = (type *)im1->image[y];  \
            for (x = 0; x < out->xsize; x++) { \
                *p0 = op(type, *p1);           \
                p0++;                          \
                p1++;                          \
            }                                  \
        }                                      \
    }

#define BINOP(name, op, type)                          \
    void name(Imaging out, Imaging im1, Imaging im2) { \
        int x, y;                                      \
        for (y = 0; y < out->ysize; y++) {             \
            type *p0 = (type *)out->image[y];          \
            type *p1 = (type *)im1->image[y];          \
            type *p2 = (type *)im2->image[y];          \
            for (x = 0; x < out->xsize; x++) {         \
                *p0 = op(type, *p1, *p2);              \
                p0++;                                  \
                p1++;                                  \
                p2++;                                  \
            }                                          \
        }                                              \
    }

#define NEG(type, v1) -(v1)
#define INVERT(type, v1) ~(v1)

#define ADD(type, v1, v2) (v1) + (v2)
#define SUB(type, v1, v2) (v1) - (v2)
#define MUL(type, v1, v2) (v1) * (v2)

#define MIN(type, v1, v2) ((v1) < (v2)) ? (v1) : (v2)
#define MAX(type, v1, v2) ((v1) > (v2)) ? (v1) : (v2)

#define AND(type, v1, v2) (v1) & (v2)
#define OR(type, v1, v2) (v1) | (v2)
#define XOR(type, v1, v2) (v1) ^ (v2)
#define LSHIFT(type, v1, v2) (v1) << (v2)
#define RSHIFT(type, v1, v2) (v1) >> (v2)

#define ABS_I(type, v1) abs((v1))
#define ABS_F(type, v1) fabs((v1))

/* --------------------------------------------------------------------
 * some day, we should add FPE protection mechanisms.  see pyfpe.h for
 * details.
 *
 * PyFPE_START_PROTECT("Error in foobar", return 0)
 * PyFPE_END_PROTECT(result)
 */

#define DIV_I(type, v1, v2) ((v2) != 0) ? (v1) / (v2) : 0
#define DIV_F(type, v1, v2) ((v2) != 0.0F) ? (v1) / (v2) : 0.0F

#define MOD_I(type, v1, v2) ((v2) != 0) ? (v1) % (v2) : 0
#define MOD_F(type, v1, v2) ((v2) != 0.0F) ? fmod((v1), (v2)) : 0.0F

static int
powi(int x, int y) {
    double v = pow(x, y) + 0.5;
    if (errno == EDOM) {
        return 0;
    }
    if (v < MIN_INT32) {
        v = MIN_INT32;
    } else if (v > MAX_INT32) {
        v = MAX_INT32;
    }
    return (int)v;
}

#define POW_I(type, v1, v2) powi(v1, v2)
#define POW_F(type, v1, v2) powf(v1, v2) /* FIXME: EDOM handling */

#define DIFF_I(type, v1, v2) abs((v1) - (v2))
#define DIFF_F(type, v1, v2) fabs((v1) - (v2))

#define EQ(type, v1, v2) (v1) == (v2)
#define NE(type, v1, v2) (v1) != (v2)
#define LT(type, v1, v2) (v1) < (v2)
#define LE(type, v1, v2) (v1) <= (v2)
#define GT(type, v1, v2) (v1) > (v2)
#define GE(type, v1, v2) (v1) >= (v2)

UNOP(abs_I, ABS_I, INT32)
UNOP(neg_I, NEG, INT32)

BINOP(add_I, ADD, INT32)
BINOP(sub_I, SUB, INT32)
BINOP(mul_I, MUL, INT32)
BINOP(div_I, DIV_I, INT32)
BINOP(mod_I, MOD_I, INT32)
BINOP(pow_I, POW_I, INT32)
BINOP(diff_I, DIFF_I, INT32)

UNOP(invert_I, INVERT, INT32)
BINOP(and_I, AND, INT32)
BINOP(or_I, OR, INT32)
BINOP(xor_I, XOR, INT32)
BINOP(lshift_I, LSHIFT, INT32)
BINOP(rshift_I, RSHIFT, INT32)

BINOP(min_I, MIN, INT32)
BINOP(max_I, MAX, INT32)

BINOP(eq_I, EQ, INT32)
BINOP(ne_I, NE, INT32)
BINOP(lt_I, LT, INT32)
BINOP(le_I, LE, INT32)
BINOP(gt_I, GT, INT32)
BINOP(ge_I, GE, INT32)

UNOP(abs_F, ABS_F, FLOAT32)
UNOP(neg_F, NEG, FLOAT32)

BINOP(add_F, ADD, FLOAT32)
BINOP(sub_F, SUB, FLOAT32)
BINOP(mul_F, MUL, FLOAT32)
BINOP(div_F, DIV_F, FLOAT32)
BINOP(mod_F, MOD_F, FLOAT32)
BINOP(pow_F, POW_F, FLOAT32)
BINOP(diff_F, DIFF_F, FLOAT32)

BINOP(min_F, MIN, FLOAT32)
BINOP(max_F, MAX, FLOAT32)

BINOP(eq_F, EQ, FLOAT32)
BINOP(ne_F, NE, FLOAT32)
BINOP(lt_F, LT, FLOAT32)
BINOP(le_F, LE, FLOAT32)
BINOP(gt_F, GT, FLOAT32)
BINOP(ge_F, GE, FLOAT32)

HPyDef_METH(unop, "unop", unop_impl, HPyFunc_VARARGS)
static HPy unop_impl(HPyContext *ctx, HPy self, HPy *args, HPy_ssize_t nargs) {
    Imaging out;
    Imaging im1;
    void (*unop)(Imaging, Imaging);

    HPy_ssize_t op, i0, i1;
    if (!HPyArg_Parse(ctx, NULL, args, nargs, "nnn", &op, &i0, &i1)) {
        return HPy_NULL;
    }

    out = (Imaging)i0;
    im1 = (Imaging)i1;

    unop = (void *)op;

    unop(out, im1);

    return HPy_Dup(ctx, ctx->h_None);
}

HPyDef_METH(binop, "binop", binop_impl, HPyFunc_VARARGS)
static HPy binop_impl(HPyContext *ctx, HPy self, HPy *args, HPy_ssize_t nargs) {
    Imaging out;
    Imaging im1;
    Imaging im2;
    void (*binop)(Imaging, Imaging, Imaging);

    Py_ssize_t op, i0, i1, i2;
    if (!HPyArg_Parse(ctx, NULL, args, nargs, "nnnn", &op, &i0, &i1, &i2)) {
        return HPy_NULL;
    }

    out = (Imaging)i0;
    im1 = (Imaging)i1;
    im2 = (Imaging)i2;

    binop = (void *)op;

    binop(out, im1, im2);

    return HPy_Dup(ctx, ctx->h_None);
}

static HPyDef *module_defines[] = {
    &unop,
    &binop,
    NULL,
};

static void
install(HPyContext *ctx, HPy h_dict, char *name, void *value) {
    HPy h_key = HPyUnicode_FromString(ctx, name);
    HPy h_value = HPyLong_FromSsize_t(ctx, (HPy_ssize_t)value);

    if (HPy_IsNull(h_key) || HPy_IsNull(h_value) || HPy_SetItem(ctx, h_dict, h_key, h_value)) {
        HPyErr_Clear(ctx);
    }
}

static int
setup_module(HPyContext *ctx, HPy h_module) {
    // TODO(hpy): This is a bit hacky
    HPy h_dict = HPy_GetAttr(ctx, h_module, HPyUnicode_FromString(ctx, "__dict__"));

    install(ctx, h_dict, "abs_I", abs_I);
    install(ctx, h_dict, "neg_I", neg_I);
    install(ctx, h_dict, "add_I", add_I);
    install(ctx, h_dict, "sub_I", sub_I);
    install(ctx, h_dict, "diff_I", diff_I);
    install(ctx, h_dict, "mul_I", mul_I);
    install(ctx, h_dict, "div_I", div_I);
    install(ctx, h_dict, "mod_I", mod_I);
    install(ctx, h_dict, "min_I", min_I);
    install(ctx, h_dict, "max_I", max_I);
    install(ctx, h_dict, "pow_I", pow_I);

    install(ctx, h_dict, "invert_I", invert_I);
    install(ctx, h_dict, "and_I", and_I);
    install(ctx, h_dict, "or_I", or_I);
    install(ctx, h_dict, "xor_I", xor_I);
    install(ctx, h_dict, "lshift_I", lshift_I);
    install(ctx, h_dict, "rshift_I", rshift_I);

    install(ctx, h_dict, "eq_I", eq_I);
    install(ctx, h_dict, "ne_I", ne_I);
    install(ctx, h_dict, "lt_I", lt_I);
    install(ctx, h_dict, "le_I", le_I);
    install(ctx, h_dict, "gt_I", gt_I);
    install(ctx, h_dict, "ge_I", ge_I);

    install(ctx, h_dict, "abs_F", abs_F);
    install(ctx, h_dict, "neg_F", neg_F);
    install(ctx, h_dict, "add_F", add_F);
    install(ctx, h_dict, "sub_F", sub_F);
    install(ctx, h_dict, "diff_F", diff_F);
    install(ctx, h_dict, "mul_F", mul_F);
    install(ctx, h_dict, "div_F", div_F);
    install(ctx, h_dict, "mod_F", mod_F);
    install(ctx, h_dict, "min_F", min_F);
    install(ctx, h_dict, "max_F", max_F);
    install(ctx, h_dict, "pow_F", pow_F);

    install(ctx, h_dict, "eq_F", eq_F);
    install(ctx, h_dict, "ne_F", ne_F);
    install(ctx, h_dict, "lt_F", lt_F);
    install(ctx, h_dict, "le_F", le_F);
    install(ctx, h_dict, "gt_F", gt_F);
    install(ctx, h_dict, "ge_F", ge_F);

    return 0;
}


HPy_MODINIT(_imagingmath)
static HPy init__imagingmath_impl(HPyContext *ctx) {

    static HPyModuleDef module_def = {
        HPyModuleDef_HEAD_INIT,
        .m_name = "_imagingmath",
        .m_doc = NULL,
        .m_size = -1,
        .defines = module_defines,
    };

    HPy h_module = HPyModule_Create(ctx, &module_def);
    if (HPy_IsNull(h_module))
        return HPy_NULL;

    if (setup_module(ctx, h_module) < 0) {
        return HPy_NULL;
    }

    return h_module;
}
