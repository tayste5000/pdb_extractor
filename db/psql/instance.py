import psycopg2
import credentials
import numpy as np
import os
import os.path as path

class Instance ( object ):

	def __init__( self ):

		self._conn = psycopg2.connect( "dbname={} user={} password={}".format( *credentials.as_tuple() ) )

		self._cursor = self._conn.cursor()

	def setup( self ):

		f = open( path.join( path.dirname( __file__ ), "init.sql" ) )

		setup_script = f.read()

		f.close()

		self._cursor.execute( setup_script )

		self._conn.commit()

	def create_dihedral( self, *args ):

		query, parameters = create_dihedral_query( *args )

		self._cursor.execute( query, parameters )

	def create_segment( self, *args ):

		query, parameters = create_segment_query( *args )

		self._cursor.execute( query, parameters )

		return self._cursor.fetchone()[0]

	def finish( self ):

		self._conn.commit()

		self._cursor.close()
		
		self._conn.close()

def create_segment_query( *args ):

	query = "INSERT INTO segment (length) VALUES (%s) RETURNING id;"

	return query, args

def create_dihedral_query( *args ):

	query = "INSERT INTO dihedral (phi,psi,res,segid,segpos) VALUES (%s,%s,%s,%s,%s);"

	return query, args