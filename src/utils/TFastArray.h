//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File MolSim.cpp
// contains main class MolSim
//------------------------------------------------------------------------------------------------

#ifndef TFASTARRAY_HEADER_
#define TFASTARRAY_HEADER_

//all util files
#include "Base.h"

#include <cstdio>
#include <cassert>

namespace utils
{

/// Array class, using C heap to allow reallocation
template<typename T> class TFastArray
{
private:
	unsigned int			_size;
	unsigned int			_capacity;
	T						*_data;
	static const unsigned int		allocation_factor = 5;	

	inline void reallocate()
	{
		// lazy initialization
		_capacity += allocation_factor;
		_data = (T*)realloc(_data, sizeof(T) * _capacity);		
	}

public:
	TFastArray()
	{
		_capacity = 0;
		_size = 0;
		_data = NULL;
	}

	/// copy constructor
	TFastArray(const TFastArray<T>& fa)
	{
		_capacity = fa._capacity;
		_size = fa._size;
		_data = NULL;
		if(_capacity != 0)
		{
			_data = realloc(_data, sizeof(T) * _capacity);

			memcpy(_data, fa._data, sizeof(T) * _capacity);
		}
		else _data = NULL;
	}

	~TFastArray()
	{
		if(_data)free(_data);
		_data = NULL;
		_size = 0; 
		_capacity = 0;
	}

	TFastArray<T>& operator = (const TFastArray<T>& fa)
	{
		_capacity = fa._capacity;
		_size = fa._size;
		_data = (T*)malloc(sizeof(T) * _capacity);
		memcpy(_data, fa._data, sizeof(T) * _capacity);
		
		return *this;
	}

	/// add element to dynamic array
	void	push_back(const T& t)
	{
		// enough space?
		if(_size >= _capacity)reallocate();		

		_data[_size] = t;
		_size++;
		
	}

	inline unsigned int size() {return _size;}
	inline unsigned int capacity()	{return _capacity;}

	inline T& operator [] (size_t index)
	{
#ifdef DEBUG
		assert(index < _size && index >= 0);
#endif
		return _data[index];
	}

	/// returns index of next element or size if it was last element
	unsigned int	erase(size_t index)
	{
#ifdef DEBUG
		assert(index < _size && index >= 0);
#endif
		// switch
		_data[index] = _data[_size - 1];
		_size--;
		return index + 1;
	}
};

}
#endif
