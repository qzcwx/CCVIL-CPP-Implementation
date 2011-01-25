#ifndef _ARCHIVE_H
#define _ARCHIVE_H

#include <EALib/PopulationT.h>

class Archive{
protected:
	unsigned MAX_NP;
	unsigned dimension;
	PopulationT<double>* pop;

public:
	Archive(unsigned inNP, unsigned inD);
	~Archive();

	unsigned getNP();
	PopulationT<double>* getPop();
};

#endif
