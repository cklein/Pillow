from tester import *

from PIL import Image, PyAccess

import test_image_putpixel as put
import test_image_getpixel as get



try:
    import cffi
except:
    skip()

Image.USE_CFFI_ACCESS = True

def test_put():
    put.test_sanity()

def test_get():
    get.test_pixel()
    get.test_image()

def _test_get_access(im):
    """ Do we get the same thing as the old pixel access """

    """ Using private interfaces, forcing a capi access and a pyaccess for the same image """
    caccess = im.im.pixel_access(False)
    access = PyAccess.new(im, False)

    w,h = im.size
    for x in range(0,w,10):
        for y in range(0,h,10):
            assert_equal(caccess[(x,y)], access[(x,y)])

def test_get_vs_c():
    _test_get_access(lena('RGB'))
    _test_get_access(lena('RGBA'))
    _test_get_access(lena('L'))
    _test_get_access(lena('LA'))
    _test_get_access(lena('1'))
    _test_get_access(lena('P'))
    #_test_get_access(lena('PA')) # PA   -- how do I make a PA image???


def _test_set_access(im, color):
    """ Are we writing the correct bits into the image? """

    """ Using private interfaces, forcing a capi access and a pyaccess for the same image """
    caccess = im.im.pixel_access(False)
    access = PyAccess.new(im, False)

    w,h = im.size
    for x in range(0,w,10):
        for y in range(0,h,10):
            access[(x,y)] = color
            assert_equal(caccess[(x,y)], color)

def test_set_vs_c():
    _test_set_access(lena('RGB'), (255, 128,0) )
    _test_set_access(lena('RGBA'), (255, 192, 128, 0))
    _test_set_access(lena('L'), 128)
    _test_set_access(lena('LA'), (128,128))
    _test_set_access(lena('1'), 255)
    _test_set_access(lena('P') , 128)
    ##_test_set_access(i, (128,128)) #PA  -- undone how to make

    
