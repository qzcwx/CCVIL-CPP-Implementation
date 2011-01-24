#ifndef _ARCHIVE_H
#define _ARCHIVE_H

class Archive{
public:
	unsigned NP;
	//TODO: Adjust the population, 
	double* pop;
	double* fitVal;

	Archive();
	~Archive();
};

#endif
