#include "Header.h"

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

	PopulationT<double> parents(runParam.NP,  ChromosomeT< double >(runParam.dimension));

	for (funcIndex = 0; funcIndex < runParam.functionToRun.size(); funcIndex++ ){
		fp = generateFuncObj(runParam, funcToRun[funcIndex]);
	}

	printf("function value = %1.20E\n", fp->compute(X));

	delete fp;
	delete[] X;

	return 0;
}


Benchmarks* generateFuncObj(RunParameter runParam, int funcID){
	Benchmarks *fp;
	// run each of specified function in "configure.ini"
	if (funcID==1){
		fp = new F1(runParam);
	}else if (funcID==2){
		fp = new F2(runParam);
	}else if (funcID==3){
		fp = new F3(runParam);
	}else if (funcID==4){
		fp = new F4(runParam);
	}else if (funcID==5){
		fp = new F5(runParam);
	}else if (funcID==6){
		fp = new F6(runParam);
	}else if (funcID==7){
		fp = new F7(runParam);
	}else if (funcID==8){
		fp = new F8(runParam);
	}else if (funcID==9){
		fp = new F9(runParam);
	}else if (funcID==10){
		fp = new F10(runParam);
	}else if (funcID==11){
		fp = new F11(runParam);
	}else if (funcID==12){
		fp = new F12(runParam);
	}else if (funcID==13){
		fp = new F13(runParam);
	}else if (funcID==14){
		fp = new F14(runParam);
	}else if (funcID==15){
		fp = new F15(runParam);
	}else if (funcID==16){
		fp = new F16(runParam);
	}else if (funcID==17){
		fp = new F17(runParam);
	}else if (funcID==18){
		fp = new F18(runParam);
	}else if (funcID==19){
		fp = new F19(runParam);
	}else if (funcID==20){
		fp = new F20(runParam);
	}else{
		cerr<<"Fail to locate Specified Function Index"<<endl;
		exit(-1);
	}
	return fp;
}
