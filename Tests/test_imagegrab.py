from helper import unittest, PillowTestCase

import sys

try:
    from PIL import ImageGrab

    class TestImageGrab(PillowTestCase):

        def test_grab(self):
            im = ImageGrab.grab()
            self.assert_image(im, im.mode, im.size)

except ImportError:
    class TestImageGrab(PillowTestCase):
        def test_skip(self):
            self.skipTest("ImportError")


class TestImageGrabImport(PillowTestCase):

    def test_import(self):
        # Arrange
        exception = None

        # Act
        try:
            from PIL import ImageGrab
            ImageGrab.__name__  # dummy to prevent Pyflakes warning
        except Exception as e:
            exception = e

        # Assert
        if sys.platform in ["win32", "darwin"]:
            self.assertIsNone(exception)
        else:
            self.assertIsInstance(exception, ImportError)
            self.assertEqual(str(exception),
                             "ImageGrab is macOS and Windows only")


if __name__ == '__main__':
    unittest.main()
