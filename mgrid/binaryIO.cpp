#include "binaryIO.h"
#include <cassert>

using namespace std;

void read_binMatrix(const string& file, vector<int> &cnt, vector<int> &col, vector<double> &ele)
{

    ifstream ifs(file, ios_base::in | ios_base::binary);

    if(!(ifs.is_open() && ifs.good()))
    {
        cerr << "ReadBinMatrix: Error cannot open file " << file << endl;
        assert(ifs.is_open());
    }
    cout << "ReadBinMatrix: Opened file " << file << endl;

    int _size;

    ifs.read(reinterpret_cast<char*>(&_size), sizeof(int));   // old: ifs.read((char*)&_size, sizeof(int));
    cnt.resize(_size);
    cout << "ReadBinMatrix: cnt size: " << _size << endl;

    ifs.read(reinterpret_cast<char*>(&_size), sizeof(int));
    col.resize(_size);
    cout << "ReadBinMatrix: col size: " << _size << endl;

    ifs.read(reinterpret_cast<char*>(&_size), sizeof(int));
    ele.resize(_size);
    cout << "ReadBinMatrix: ele size: " << _size << endl;


    ifs.read(reinterpret_cast<char*>(cnt.data()), cnt.size() * sizeof(int));
    ifs.read(reinterpret_cast<char*>(col.data()), col.size() * sizeof(int));
    ifs.read(reinterpret_cast<char*>(ele.data()), ele.size() * sizeof(double));

    ifs.close();
    cout << "ReadBinMatrix: Finished reading matrix.." << endl;

}

void write_binMatrix(const string& file, const vector<int> &cnt, const vector<int> &col, const vector<double> &ele)
{
    ofstream ofs(file, ios_base::out | ios_base::binary);


    if(!(ofs.is_open() && ofs.good()))
    {
        cerr << "WriteBinMatrix: Error cannot open file " << file << endl;
        assert(ofs.is_open());
    }
    cout << "WriteBinMatrix: Opened file " << file << endl;

    int _size = static_cast<int>( cnt.size() );
    cout << "WriteBinMatrix: cnt size: " << _size << endl;
    ofs.write(reinterpret_cast<char*>(&_size), sizeof(int));
    _size = static_cast<int>( col.size() );
    cout << "WriteBinMatrix: col size: " << _size << endl;
    ofs.write(reinterpret_cast<char*>(&_size), sizeof(int));
    _size = static_cast<int>( ele.size() );
    cout << "WriteBinMatrix: ele size: " << _size << endl;
    ofs.write(reinterpret_cast<char*>(&_size), sizeof(int));

    ofs.write(reinterpret_cast<char const *>(cnt.data()), cnt.size() * sizeof(int));
    ofs.write(reinterpret_cast<char const *>(col.data()), col.size() * sizeof(int));
    ofs.write(reinterpret_cast<char const *>(ele.data()), ele.size() * sizeof(double));

    ofs.close();
    cout << "WriteBinMatrix: Finished writing matrix.." << endl;
}
