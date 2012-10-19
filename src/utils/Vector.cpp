// added by studentls

#include "Vector.h"


template <typename type, int length>
std::ostream& operator << (std::ostream& stream, const utils::Vector<type, length>& v) {

	stream << "[";
	for (int i = 0; i < length; i++) {
		stream << v.content[i] << ";";
	}
	stream << "]";
	return stream;
}