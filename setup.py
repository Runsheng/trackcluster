from setuptools import setup, Extension

setup(
    name='trackcluster',
    version='0.1.0',
    packages=['', 'ssw', 'test', 'script'],
    ext_modules= [Extension("ssw.trackcluster", ["ssw/ssw.c"])],
    url='https://github.com/runsheng/trackcluster',
    license='GPL-2',
    author='runsheng',
    author_email='runsheng.lee@gmail.com',
    description='RNA-seq analysis for Nanopore direct-RNA sequencing'
)
