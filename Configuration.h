//===========================================================================
//===========================================================================


#include <FileUtil/Params.h>


class MyParams : public Params
{
public:
    MyParams(int argc, char **argv) : Params(argc, argv)
    {
        setDefault();
        scanFrom(confFile);
    }

    void io(std::istream& is, std::ostream& os, FileUtil::iotype type)
    {
        FileUtil::io_strict(is, os, "n"        , n     , 4u     , type);
        FileUtil::io_strict(is, os, "d"        , d     , 3.14   , type);
        FileUtil::io_strict(is, os, "f"        , f     , 3.14   , type);
    }

    // Method to show the current content of the class variable:
    void monitor(std::ostream &os) const
    {
		os << "n\t"       << n        << "\t(the n)\n"
			<< "d\t"       << d        << "\t(the d)\n"
			<< "f\t"       << f        << "\t(the f)\n" << std::endl;
    }

    friend std::ostream& operator <<(std::ostream &os, const MyParams &obj);

    unsigned n;
    double   d;
    double   f;
};

std::ostream& operator <<(std::ostream &os, const MyParams &obj)
{
    obj.monitor(os);
    return os;
}

//
// call ./fileUtilClass -conf parameter.conf
//
int main(int argc, char **argv)
{
	MyParams p(argc, argv);
	std::cout << p << std::endl;
}
