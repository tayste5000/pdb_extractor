DROP TABLE IF EXISTS segment CASCADE;
DROP TABLE IF EXISTS dihedral;

CREATE TABLE segment(
	id SERIAL,
	length INT NOT NULL,
	PRIMARY KEY(id)
);

CREATE TABLE dihedral(
	id SERIAL,
	phi double precision NOT NULL,
	psi double precision NOT NULL,
	res CHAR(1) NOT NULL,
	segid INT NOT NULL,
	segpos INT NOT NULL,
	PRIMARY KEY(id),
	FOREIGN KEY(segid) REFERENCES segment (id) ON DELETE CASCADE
);