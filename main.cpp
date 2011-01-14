#include "Header.h"
#include "RunParameter.h"

#include "F1.h"
#include "F2.h"
#include "F3.h"
#include "F4.h"
#include "F5.h"
#include "F6.h"
#include "F7.h"
#include "F8.h"
#include "F9.h"
#include "F10.h"
#include "F11.h"
#include "F12.h"
#include "F13.h"
#include "F14.h"
#include "F15.h"
#include "F16.h"
#include "F17.h"
#include "F18.h"
#include "F19.h"
#include "F20.h"

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
		}else if (funcToRun[funcIndex]==5){
			fp = new F5(runParam);
		}else if (funcToRun[funcIndex]==6){
			fp = new F6(runParam);
		}else if (funcToRun[funcIndex]==7){
			fp = new F7(runParam);
		}else if (funcToRun[funcIndex]==8){
			fp = new F8(runParam);
		}else if (funcToRun[funcIndex]==9){
			fp = new F9(runParam);
		}else if (funcToRun[funcIndex]==10){
			fp = new F10(runParam);
		}else if (funcToRun[funcIndex]==11){
			fp = new F11(runParam);
		}else if (funcToRun[funcIndex]==12){
			fp = new F12(runParam);
		}else if (funcToRun[funcIndex]==13){
			fp = new F13(runParam);
		}else if (funcToRun[funcIndex]==14){
			fp = new F14(runParam);
		}else if (funcToRun[funcIndex]==15){
			fp = new F15(runParam);
		}else if (funcToRun[funcIndex]==16){
			fp = new F16(runParam);
		}else if (funcToRun[funcIndex]==17){
			fp = new F17(runParam);
		}else if (funcToRun[funcIndex]==18){
			fp = new F18(runParam);
		}else if (funcToRun[funcIndex]==19){
			fp = new F19(runParam);
		}else if (funcToRun[funcIndex]==20){
			fp = new F20(runParam);
		}else{
			cerr<<"Fail to locate Specified Function Index"<<endl;
			exit(-1);
		}
	}

	printf("function value = %1.20E\n", fp->compute(X));

	delete[] X;
	return 0;
}


