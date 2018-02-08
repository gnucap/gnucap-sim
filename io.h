
inline OMSTREAM& operator<<(OMSTREAM& o, COMPLEX c)
{
	o << c.real();
	if(c.imag() <0){
		o << "-" << -c.imag();
	}else{
		o << "+" << c.imag();
	}
	return  o	<< "* i";
}
