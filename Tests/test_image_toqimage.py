from helper import unittest, PillowTestCase, hopper
from test_imageqt import PillowQtTestCase

from PIL import ImageQt


if ImageQt.qt_is_installed:
    from PIL.ImageQt import QImage


class TestToQImage(PillowQtTestCase, PillowTestCase):

    def test_sanity(self):
        PillowQtTestCase.setUp(self)
        for mode in ('1', 'RGB', 'RGBA', 'L', 'P'):
            data = ImageQt.toqimage(hopper(mode))

            self.assertTrue(isinstance(data, QImage))
            self.assertFalse(data.isNull())

            # Test saving the file
            tempfile = self.tempfile('temp_{}.png'.format(mode))
            data.save(tempfile)


if __name__ == '__main__':
    unittest.main()

# End of file
