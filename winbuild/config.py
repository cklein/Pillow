import os

SF_MIRROR = 'http://iweb.dl.sourceforge.net'

pythons = {  # '26': 7,
           '27': 7,
           # '32': 7,
           '33': 7.1,
           '34': 7.1}

VIRT_BASE = "c:/vp/"
X64_EXT = os.environ.get('X64_EXT', "x64")

libs = {
    'zlib': {
        'url': 'http://zlib.net/zlib128.zip',
        'hash': 'md5:126f8676442ffbd97884eb4d6f32afb4',
        'dir': 'zlib-1.2.8',
    },
    'jpeg': {
        'url': 'http://www.ijg.org/files/jpegsr9b.zip',
        'hash': 'md5:a21b8024d78ba05857a75272b4fa95ec',  # not found - generated by wiredfool
        'dir': 'jpeg-9b',
    },
    'tiff': {
        'url': 'ftp://ftp.remotesensing.org/pub/libtiff/tiff-4.0.6.zip',
        'hash': 'md5:f5b485d750b2001255ed64224b98b857',
        'dir': 'tiff-4.0.6',
    },
    'freetype': {
        'url': 'http://download.savannah.gnu.org/releases/freetype/freetype-2.6.3.tar.gz',
        'hash': 'md5:8b0c80b64042b7c39b9fa9debe6305c3',
        'dir': 'freetype-2.6.3',
    },
    'lcms': {
        'url': SF_MIRROR+'/project/lcms/lcms/2.7/lcms2-2.7.zip',
        'hash': 'sha1:7ff1a5b721ca719760ba6eb4ec6f38d5e65381cf',
        'dir': 'lcms2-2.7',
    },
    'tcl-8.5': {
        'url': SF_MIRROR+'/project/tcl/Tcl/8.5.18/tcl8518-src.zip',
        'hash': 'sha1:4c2aed9043088c630a4c795265e2738ef1b4db3b',
        'dir': '',
    },
    'tk-8.5': {
        'url': SF_MIRROR+'/project/tcl/Tcl/8.5.18/tk8518-src.zip',
        'hash': 'sha1:273f55148777413774aa722ecad25cabda1e31ae',
        'dir': '',
        'version': '8.5.18',
    },
    'tcl-8.6': {
        'url': SF_MIRROR+'/project/tcl/Tcl/8.6.4/tcl864-src.zip',
        'hash': 'md5:35748d2fc61e08a2fdb23b85c6f8c4a0',
        'dir': '',
    },
    'tk-8.6': {
        'url': SF_MIRROR+'/project/tcl/Tcl/8.6.4/tk864-src.zip',
        'hash': 'md5:111d45061a69e7f5250b6ec8ca7c4f35',
        'dir': '',
        'version': '8.6.4',
    },
    'webp': {
        'url': 'http://downloads.webmproject.org/releases/webp/libwebp-0.5.0.tar.gz',
        'hash': 'sha1:d3de815b272fcf88fc4f2dc1ab65d176bcb8df68',
        'dir': 'libwebp-0.5.0',
    },
    'openjpeg': {
        'url': SF_MIRROR+'/project/openjpeg/openjpeg/2.1.0/openjpeg-2.1.0.tar.gz',
        'hash': 'md5:f6419fcc233df84f9a81eb36633c6db6',
        'dir': 'openjpeg-2.1.0',
    },
}

bin_libs = {
    'openjpeg': {
        'filename': 'openjpeg-2.0.0-win32-x86.zip',
        'hash': 'sha1:xxx',
        'version': '2.0'
    },
}

compilers = {
    (7, 64): {
        'env_version': 'v7.0',
        'vc_version': '2008',
        'env_flags': '/x64 /xp',
        'inc_dir': 'msvcr90-x64',
        'platform': 'x64',
        'webp_platform': 'x64',
    },
    (7, 32): {
        'env_version': 'v7.0',
        'vc_version': '2008',
        'env_flags': '/x86 /xp',
        'inc_dir': 'msvcr90-x32',
        'platform': 'Win32',
        'webp_platform': 'x86',
    },
    (7.1, 64): {
        'env_version': 'v7.1',
        'vc_version': '2010',
        'env_flags': '/x64 /vista',
        'inc_dir': 'msvcr10-x64',
        'platform': 'x64',
        'webp_platform': 'x64',
    },
    (7.1, 32): {
        'env_version': 'v7.1',
        'vc_version': '2010',
        'env_flags': '/x86 /vista',
        'inc_dir': 'msvcr10-x32',
        'platform': 'Win32',
        'webp_platform': 'x86',
    },
}


def pyversion_fromEnv():
    py = os.environ['PYTHON']

    py_version = '27'
    for k in pythons.keys():
        if k in py:
            py_version = k
            break

    if '64' in py:
        py_version = '%s%s' % (py_version, X64_EXT)

    return py_version


def compiler_fromEnv():
    py = os.environ['PYTHON']

    for k, v in pythons.items():
        if k in py:
            compiler_version = v
            break

    bit = 32
    if '64' in py:
        bit = 64

    return compilers[(compiler_version, bit)]
