# phipsidb
Simple python script for getting dihedral angles from PDB structures and storing them in a database

### Prerequisites:

Install correct python driver for your database (e.g. for postgresql install psycopg2)

Create a database before running

In the db/ directory, find the folder corresponding to the database you are using, copy the contents of the credentials.example.py file to credentials.py, and replace the placholders with actual database login info

### Database support

postgresql - psql

The code for interfacing with the database is fairly simple so feel free to submit a pull request adding support for the database of your choice or submit an issue requesting support and I should be able to add it

### How to use:

    cd phipsidb
    python process_pdb.py databasetype numberofpdbfiles

for example
  
    python process_pdb.py psql 10000
    
will pull 10,000 PDB structures, compute all of the dihedrals angles, and store them in the postgresql database you have specified in db/psql/credentials.py

The script will handle creating database tables to store your data in, however, be careful if you are using a database that already has data, the script will drop any records in tables named "dihedral" or "segment" prior to running
