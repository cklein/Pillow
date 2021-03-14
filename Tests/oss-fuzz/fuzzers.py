import warnings
import io

from PIL import Image, ImageFont, ImageDraw, ImageFilter, ImageFile, PcfFontFile

def enable_decompressionbomb_error():
    ImageFile.LOAD_TRUNCATED_IMAGES = True
    warnings.filterwarnings("ignore")
    warnings.simplefilter("error", Image.DecompressionBombWarning)

def fuzz_image(data):
    # This will fail on some images in the corpus, as we have many
    # invalid images in the test suite.
    with Image.open(io.BytesIO(data)) as im:
        im.rotate(45)
        im.filter(ImageFilter.DETAIL)
        im.save(io.BytesIO(), "BMP")

def fuzz_font(data):
    # This should not fail on a valid font load for any of the fonts in the corpus
    wrapper = io.BytesIO(data)
    try:
        font = ImageFont.truetype(wrapper)
    except OSError:
        # pcf/pilfonts/random garbage here here. They're different.
        return

    font.getsize_multiline("ABC\nAaaa")
    font.getmask("test text")
    with Image.new(mode="RGBA", size=(200, 200)) as im:
        draw = ImageDraw.Draw(im)
        draw.multiline_textsize("ABC\nAaaa", font, stroke_width=2)
        draw.text((10,10), "Test Text", font=font, fill="#000")
