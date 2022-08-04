from setuptools import setup, Extension, find_packages
from trackcluster import __version__

setup(
    name='trackcluster',
    version=__version__,
    packages=['', "trackcluster", 'test', 'script'],
    #ext_modules= [Extension("ssw.trackcluster", ["ssw/ssw.c"])],
    url='https://github.com/runsheng/trackcluster',
    license='GPL-2',
    author='runsheng',
    author_email='runsheng.lee@gmail.com',
    description='RNA-seq analysis for Long read RNA sequencing',

    install_requires=["biopython>=1.78",
                      "numpy>=1.19.4",
                      "pandas>=1.1.5",
                      "pysam>=0.16.0.1",
                      "tqdm>=4.64.0"],

    scripts=['script/trackrun.py',
             'script/bam2bigg.py',
             'script/bigg2b.py',
             'script/gff2bigg.py',
             'script/biggmutant.py'],

)
