#include <iostream>

template <class T>
vector<T> vcat(const vector<T>& top, const vector<T>& bottom)
{
    vector<T> ret = top;

    // Option 1: insert()
    ret.insert(ret.end(), bottom.begin(), bottom.end());

    return ret;
}

int main(){
	

}