#include "Header.h"
#include "RunParameter.h"
#include "F1.h"
#include "F2.h"
#include "F3.h"
#include "F4.h"

int main(){
	RunParameter runParam;
	double* X;
	unsigned int funcIndex;
	Benchmarks *fp;
	vector<int> funcToRun=runParam.functionToRun;

	X = new double[runParam.dimension];

	for (int i=0; i<runParam.dimension; i++){
		X[i]=0;
	}

	for (funcIndex = 0; funcIndex < runParam.functionToRun.size(); funcIndex++ ){
		// run each of specified function in "configure.ini"
		if (funcToRun[funcIndex]==1){
			fp = new F1(runParam);
		}else if (funcToRun[funcIndex]==2){
			fp = new F2(runParam);
		}else if (funcToRun[funcIndex]==3){
			fp = new F3(runParam);
		}else if (funcToRun[funcIndex]==4){
			fp = new F4(runParam);
		}else{
			cerr<<"Fail to locate Specified Function Index"<<endl;
		}
	}

	printf("%1.20E\n", fp->compute(X));

	delete[] X;
	return 0;

}


