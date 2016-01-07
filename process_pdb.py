from Bio import PDB
from Bio import SeqUtils

import urllib2
import db
import os
import gzip
from StringIO import StringIO
import time
import xml.etree.ElementTree as ET
import sys
import importlib

'''
Database adapter specified by 1st command line arguments

available options:
-psql (postgresql)
'''

try:
	db = importlib.import_module( "db.{}".format( sys.argv[1] ) )

except ImportError:
	sys.exit( "Database module {} not found".format( sys.argv[1] ) )

db_instance = db.Instance()

'''
testing if a residue is a valid amino acid
'''

aas = ( "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y" )

def is_aa( aa ):

	if aa in aas:

		return True

	else:

		return False


class BackboneError( Exception ):

	'''
	Throw when there is an error parsing dihedral angles
	'''

	def __init__( self, message, res_number ):
		self.message = message
		self.res_number = res_number

	def __str__( self ):
		return repr( "BackboneError occured at residue {number}: {message}".format(
			number = self.res_number,
			message = self.message
		))


class NoHetAtm( PDB.StructureBuilder.StructureBuilder ):

	"""
	Alternative structure builder that rejects all non-peptide residues
	"""

	def init_res( self, resname, field, resseq, icode ):

		if field != " ":

			raise PDBConstructionException( "Rejecting all non-peptide residues" )

		super( NoHetAtm, self ).init_res( self, resname, field, resseq, icode )	


def check_dict( test_dict ):

	return False not in test_dict.values()


def get_dihedral( residue_list ):

	'''
	returns phi and psi angles of a residue and the amino acid sidechain present

	residue_list - []Bio.PDB.Residue - list of 3 *hopefully* continuous residues

	'''

	for one, two in zip( residue_list[:-1], residue_list[1:] ):

		if ( two.get_id()[1] - one.get_id()[1] ) != 1:

			raise BackboneError( "Discontinuous residues", two.get_id()[1] )

	atoms = (
		{"C": False},
		{"N": False,
		"CA": False,
		"C": False},
		{"N": False}
	)

	for i, residue in enumerate( residue_list ):

		if i == 1:

			res_name = SeqUtils.seq1( residue.get_resname() )

			if not is_aa( res_name ):

				raise BackboneError( "Not a valid amino acid", residue.get_id()[1] )

		for atom in residue.get_unpacked_list():

			if atom.name in atoms[i].keys():
				
				atoms[i][ atom.name ] = atom.get_vector()

	if False in map( check_dict, atoms ):

		raise BackboneError( "Missing backbone atoms", residue.get_id()[1] )

	dihedrals = [
		PDB.calc_dihedral( atoms[0]["C"], atoms[1]["N"], atoms[1]["CA"], atoms[1]["C"] ), #phi
		PDB.calc_dihedral( atoms[1]["N"], atoms[1]["CA"], atoms[1]["C"], atoms[2]["N"] ) #psi
	]

	return ( dihedrals, res_name )


def store_dihedrals( dihedrals, res_names, length ):

	segid = db_instance.create_segment( length )

	for i, (dihedral, res_name) in enumerate( zip( dihedrals, res_names ) ):

		db_instance.create_dihedral(
			dihedral[0],
			dihedral[1],
			res_name,
			segid,
			i
		)


def extract_file( code ):

	'''
	download PDB file by ID,
	extract into Bio.PDB object,
	parse dihedrals and store in database
	'''

	url = "ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/{short}/pdb{code}.ent.gz".format(
		short = code[1:3].lower() ,
		code = code.lower()
	)

	gzipped = urllib2.urlopen( url )

	buf = StringIO( gzipped.read() )

	pdb_file = gzip.GzipFile( fileobj = buf )

	gzipped.close()

	builder = NoHetAtm()

	parser = PDB.PDBParser( structure_builder = builder )

	structure = parser.get_structure( code , pdb_file )

	pdb_file.close()

	for chain in structure.get_chains():

		dihedrals = []
		res_names = []

		start = 0

		try:

			last_res = chain.get_residues().next()
			this_res = chain.get_residues().next()
			next_res = chain.get_residues().next()

		except StopIteration:

			continue

		try:

			dihedral, res_name = get_dihedral([
				last_res,
				this_res,
				next_res ])

			dihedrals.append( dihedral )
			res_names.append( res_name )

		except BackboneError:

			start = 1

		for i, res in enumerate( chain.get_residues() ):

			try:

				dihedral, res_name = get_dihedral([
					last_res,
					this_res,
					next_res ])

				dihedrals.append( dihedral )
				res_names.append( res_name )

			except BackboneError:

				length = ( i ) - start

				start = ( i + 1 )

				if length > 0:

					store_dihedrals(
						dihedrals,
						res_names,
						length
					)

					dihedrals = []
					res_names = []

			last_res, this_res, next_res = this_res, next_res, res

		length = ( i ) - start

		if length > 0:

			store_dihedrals(
				dihedrals,
				res_names,
				length
			)


def extract_many( num = 10 ):

	'''
	get list of all PDB ID's,
	select a given number from this list,
	run extract_file on each PDB
	'''

	print "\tFetching {num} pdb ID's\n".format( num = num )

	before = time.time()

	raw_xml = urllib2.urlopen( "http://www.rcsb.org/pdb/rest/getCurrent" )

	pdb_ids_xml = ET.fromstring( raw_xml.read() )

	raw_xml.close()

	pdb_ids = []

	for element in pdb_ids_xml[:num]:

		pdb_ids += element.attrib.values()

	after = time.time()

	print "\tRetrieved {num} pdb ID's in {time} s\n".format( num = num, time = after - before )

	print "\tBegining dihedral extraction ...\n"

	for code in pdb_ids:

		print "\tExtracting dihedrals from {code} ...\n".format( code = code )

		before = time.time()

		extract_file( code )

		after = time.time()

		print "\tExtracted dihedrals in {time} s\n".format( time = after - before )



db_instance.setup()

extract_many( num = int( sys.argv[2] ) )

db_instance.finish()