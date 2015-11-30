from distutils.core import setup
import glob

version = '1.2.0'
doc_prefix = 'doc/NZMATH-%s/' % version
data_prefix = 'nzmath'

setup(
    name='NZMATH',
    version=version,
    url='http://tnt.math.se.tmu.ac.jp/nzmath/',
    author='NZMATH development group',
    author_email='nzmath-user@tnt.math.se.tmu.ac.jp',
    description='number theory oriented calculation system',
    classifiers=['Development Status :: 5 - Production/Stable',
                 'License :: OSI Approved :: BSD License',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python',
                 'Topic :: Scientific/Engineering :: Mathematics',
                ],

    packages=['nzmath', 'nzmath.config', 'nzmath.factor', 'nzmath.poly', 'nzmath.plugin', 'nzmath.plugin.math'],

    data_files=[
        (data_prefix, ['data/discriminant.csv']),
        (doc_prefix + 'manual',
            glob.glob('manual/*.pdf')),
    ]
)
