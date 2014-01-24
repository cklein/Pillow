from tester import *

from PIL import Image
import os

base = 'Tests/images/bmp/'


def get_files(d, ext='.bmp'):
    return [os.path.join(base,d,f) for f
            in os.listdir(os.path.join(base, d)) if ext in f]

def xtest_bad():
    """ These shouldn't crash/dos, but they shouldn't return anything either """
    for f in get_files('b'):
        try:
            print ("Trying %s"%f)
            im = Image.open(f)
            print ("%s, %s" %(im.size, im.mode))
            im.load()
        except Exception as msg:
            print ("Bad Image %s: %s" %(f,msg))

def xtest_questionable():
    """ These shouldn't crash/dos, but its not well defined that these are in spec """
    for f in get_files('q'):
        try:
            print ("Trying %s"%f)
            im = Image.open(f)
            print ("%s, %s" %(im.size, im.mode))
            im.load()
        except Exception as msg:
            print ("Bad Image %s: %s" %(f,msg))


def test_good():
    """ These should all work. There's a set of target files in the
    html directory that we can compare against. """

    # Target files, if they're not just replacing the extension
    file_map = { 'pal1wb.bmp': 'pal1.png',
                 'pal4rle.bmp': 'pal4.png', 
                 'pal8-0.bmp': 'pal8.png',
                 'pal8rle.bmp': 'pal8.png',
                 'pal8topdown.bmp': 'pal8.png',
                 'pal8nonsquare.bmp': 'pal8nonsquare-v.png',
                 'pal8os2.bmp': 'pal8.png',
                 'pal8os2sp.bmp': 'pal8.png',
                 'pal8os2v2.bmp': 'pal8.png',
                 'pal8os2v2-16.bmp': 'pal8.png',
                 'pal8v4.bmp': 'pal8.png',
                 'pal8v5.bmp': 'pal8.png',
                 'rgb16-565pal.bmp': 'rgb16-565.png',
                 'rgb24pal.bmp': 'rgb24.png',
                 'rgb32.bmp': 'rgb24.png',
                 'rgb32bf.bmp': 'rgb24.png'
                 }     
                 
    def get_compare(f):
        (head, name) = os.path.split(f)
        if name in file_map:
            return os.path.join(base, 'html', file_map[name])
        (name,ext) = os.path.splitext(name)
        return os.path.join(base, 'html', "%s.png"%name)
    
    for f in get_files('g'):
        try:
            print '.'
            #print ("Trying %s"%f)
            im = Image.open(f)
            #print ("%s, %s" %(im.size, im.mode))
            im.load()
            compare = Image.open(get_compare(f))
            compare.load()
            im = im.convert(compare.mode)
            #print ("%s, %s" %(compare.size, compare.mode))
            if not assert_image_similar(im, compare,10):
                print (f)
                im.show()
                compare.show()
                sys.exit()
            
            
        except Exception as msg:
            print ("Bad Image %s: %s" %(f,msg))
            pass
