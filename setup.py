from setuptools import setup

setup(name='catalog_org',
	version='0.1.0',
   	author='Paula Ferreira',
   	author_email='psilf12@gmail.com',
   	packages=['bao_pz',],
   	#scripts=['bin/script1','bin/script2'],
   	#url='http://pypi.python.org/pypi/PackageName/',
   	license='LICENSE.txt',
   	description='An awesome package that does something',
  	#long_description=open('README.txt').read(),
  	install_requires=[
	"astropy",
	"numpy",
	'scipy',
	'pandas',
	'Corrfunc',],)
