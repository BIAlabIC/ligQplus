from distutils.core import setup
setup(name='uniprot_tools',
		version='0.4.6',
		py_modules=['uniprot'],
		scripts=['uniprot'],

		requires=['requests','argparse'],

		author='Jan Rudolph',
		license='MIT license',
		author_email='jan.daniel.rudolph@gmail.com',
		description='Simple interface for uniprot.org',
		url='https://github.com/jdrudolph/uniprot',
		long_description=open('README').read(),
		)
