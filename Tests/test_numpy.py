from tester import *

from PIL import Image

try:
    import site
    import numpy
except ImportError:
    skip()

def test_numpy_to_image():

    def to_image(dtype, bands=1, bool=0):
        if bands == 1:
            if bool:
                data = [0, 1] * 50
            else:
                data = list(range(100))
            a = numpy.array(data, dtype=dtype)
            a.shape = 10, 10
            i = Image.fromarray(a)
            if list(i.getdata()) != data:
                print("data mismatch for", dtype)
        else:
            data = list(range(100))
            a = numpy.array([[x]*bands for x in data], dtype=dtype)
            a.shape = 10, 10, bands
            i = Image.fromarray(a)
            if list(i.split()[0].getdata()) != list(range(100)):
                print("data mismatch for", dtype)
        # print dtype, list(i.getdata())
        return i

    # assert_image(to_image(numpy.bool, bool=1), "1", (10, 10))
    # assert_image(to_image(numpy.bool8, bool=1), "1", (10, 10))

    assert_exception(TypeError, lambda: to_image(numpy.uint))
    assert_image(to_image(numpy.uint8), "L", (10, 10))
    assert_exception(TypeError, lambda: to_image(numpy.uint16))
    assert_exception(TypeError, lambda: to_image(numpy.uint32))
    assert_exception(TypeError, lambda: to_image(numpy.uint64))

    assert_image(to_image(numpy.int8), "I", (10, 10))
    if Image._ENDIAN == '<': # Little endian
        assert_image(to_image(numpy.int16), "I;16", (10, 10))
    else:
        assert_image(to_image(numpy.int16), "I;16B", (10, 10))
    assert_image(to_image(numpy.int32), "I", (10, 10))
    assert_exception(TypeError, lambda: to_image(numpy.int64))

    assert_image(to_image(numpy.float), "F", (10, 10))
    assert_image(to_image(numpy.float32), "F", (10, 10))
    assert_image(to_image(numpy.float64), "F", (10, 10))

    assert_image(to_image(numpy.uint8, 3), "RGB", (10, 10))
    assert_image(to_image(numpy.uint8, 4), "RGBA", (10, 10))


# based on an erring example at http://is.gd/6F0esS
def test_3d_array():
    a = numpy.ones((10, 10, 10), dtype=numpy.uint8)
    assert_image(Image.fromarray(a[1, :, :]), "L", (10, 10))
    assert_image(Image.fromarray(a[:, 1, :]), "L", (10, 10))
    assert_image(Image.fromarray(a[:, :, 1]), "L", (10, 10))


def test_16bit():
    img = Image.open('Tests/images/12bit.cropped.tif')
    px = img.load()
    np_img = numpy.array(img)
    assert_equal(np_img.shape, (64,64))
    assert_equal(px[1,1],np_img[1,1])
    
