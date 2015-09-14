#
# Copyright John Reid 2007
#

"""
Draws logos of PSSMs with optional transparencies. The main entry point is L{pssm_as_image}.
"""


from PIL import Image, ImageColor, ImageFont, ImageDraw, ImageOps
import math, warnings, sys, os

_default_font = 'arial.ttf'
"The default font."

_default_font_size = 800
"The default font size."

for dir in [
    '/home/reid/.fonts',
    '/home/john/.fonts',
]:
    font_filename = os.path.join(dir, 'arial.ttf')
    if os.path.exists(font_filename):
        _default_font = font_filename




def get_font(font_name, font_size):
    """
    Tries to load the named font at the given size but falls back on default font if this is not possible.

    @return: The named font at the given size if possible.
    """
    try:
        return ImageFont.truetype(font_name, font_size)
    except:
        #print sys.exc_info()
        warnings.warn('Could not load font from %s of size %d; Trying free type font.' % (font_name, font_size))
        try:
            return ImageFont.FreeTypeFont(font_name, font_size)
        except:
            #print sys.exc_info()
            warnings.warn('Could not load free type font from %s of size %d; Trying default font.' % (font_name, font_size))
        return ImageFont.load_default()




def text_as_image(text, colour='black', font_name=_default_font, use_font_size=_default_font_size):
    """
    Creates an image of the given text with no border around the characters.

@arg text: The text.
@arg colour: The colour of the text in the image.
@arg font_name: The font name.
@arg use_font_size: The font size.
@return: A PIL image of the given text.
    """
    font = get_font(font_name, use_font_size)
    text_size = font.getsize(text)
    uncropped = Image.new("RGBA", text_size)
    draw = ImageDraw.Draw(uncropped)
    draw.text((0,0), text, font = font, fill = colour)
    cropped = uncropped.crop(uncropped.getbbox())
    result = Image.new("RGB", cropped.size, 'white')
    result.paste(cropped, None, cropped)
    return result




def entropy(p, base = 2.):
    "@return: The entropy of the distribution, p, using logarithms of the given base, b."
    try:
        if p > 0.0:
            return -p*math.log(p,base)
    except:
        print 'Math domain error: math.log(%s, %f)' % (str(p), base)
        pass
    return 0.0




def dist_as_image(dist, size):
    """
    Creates an image representing the distribution over A, C, G, and T needed for a logo of a PSSM.
    See U{http://www.lecb.ncifcrf.gov/~toms/paper/logopaper/paper/index.html}.

    @arg dist: An array length 4 of probabilities.
    @arg size: A tuple (width, height) for the size of the image.
    @return: A PIL image.
    """
    if sum(dist) - 1. > 1.e-4:
        raise ValueError('distribution sums > 1.')
    R = 2.0 - sum([ entropy(p) for p in dist ])
    information = [ (R * p, n) for p, n in zip(dist, nucleotides) ]
    information.sort(lambda e1, e2: cmp(e1[0], e2[0]))
    result = Image.new("RGB", size, 'white')
    y = size[1]
    for i, n in information:
        if i > 0.0:
            base_size = (size[0],int(size[1]*i/2))
            base_image = n.copy().resize(base_size)
            result.paste(base_image, (0,y-base_size[1]))
            y -= base_size[1]
    return result




def pssm_as_image(pssm, size=None, transparencies=None):
    """
    Creates a logo for a PSSM. Handles transparent bases.

    @arg pssm: The PSSM frequencies. I.e. a (K,4) shaped array for a PSSM with K bases.
    @arg size: The desired size of the image as a tuple of (width, height). Default will be (80*K, 240).
    @arg transparencies: An array of length K. Each entry should be between 0 and 1 and gives the transparency
    of the corresponding base. The value will be converted into a percentage and drawn on to the base's logo.
    @return: A PIL image object of the logo representing the PSSM.
    """
    if None == size:
        size = (len(pssm)*80,240)
    if None != transparencies and len(transparencies) != len(pssm):
        raise RuntimeError('transparencies array is wrong length')
    result = Image.new("RGB", size, 'white')
    base_size = (size[0]/len(pssm),size[1])

    # for each base in the pssm
    for i, dist in enumerate(pssm):
    # get the distribution at this base as an image
        dist_image = dist_as_image(dist, base_size)

        # where this base will be pasted
        location = (i*base_size[0],0)

        # do we have transparency information?
        if None != transparencies:
            mask = Image.new("L", dist_image.size, 255 * transparencies[i])

            # paste into result image
            result.paste(dist_image, box=location, mask=mask)

            # write the percentage on
            if transparencies[i] < 1.0:
                text = '%2.0f%%' % (100.0 * transparencies[i])
                font = get_font(_default_font, int(100*size[1]/800))
                draw = ImageDraw.Draw(result)
                textsize = draw.textsize(text, font)
                xy = (location[0]+base_size[0]/2-textsize[0]/2, location[1]+size[1]/2-textsize[1]/2)
                draw.text(
                        xy,
                        text=text,
                        font=font,
                        fill='black'
                )

        else:
                # no transparency info
            result.paste(dist_image, box=location)

    return result

colours = {
  'orange' : '#ffa500',
  'red'    : '#ff0000',
  'green'  : '#008000',
  'blue'   : '#0000ff',
}
A = text_as_image( 'A', colour=colours['green'])
"An image of the letter A"

C = text_as_image( 'C', colour=colours['blue'])
"An image of the letter C"

G = text_as_image( 'G', colour=colours['orange'])
"An image of the letter G"

T = text_as_image( 'T', colour=colours['red'])
"An image of the letter T"

nucleotides = [ A, C, G, T ]
"An array of images of the bases A, C, G, T"





if '__main__' == __name__:

            # for im in [A,C,G,T]: im.show()

    _dist_image = dist_as_image([.80,.01,.04,.15], (400,1200))
    #_dist_image.show()
    #dist_as_image([0,1.,0,0], (400,1200)).show()

    pssm_as_image(
            [
              [ 1,0,0,0 ],
              [ 0,.5,.5,0 ],
              [ 0,1,0,0 ],
              [ 1,0,0,0 ],
              [ 1,0,0,0 ],
            ],
            transparencies = [
              1.,
              1.,
              .5,
              .998,
              .09
            ]
    ).show()
