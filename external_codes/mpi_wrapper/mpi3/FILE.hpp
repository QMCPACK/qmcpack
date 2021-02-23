#ifndef ALF_BOOST_MPI3_FILE_HPP
#define ALF_BOOST_MPI3_FILE_HPP

namespace boost{
namespace mpi3{

struct FILE{
	MPI_File impl_;
	static void delete_(std::string const& filename){MPI_File_delete(filename.c_str(), MPI_INFO_NULL);}
	int amode() const{
		int ret;
		MPI_File_get_amode(impl_, &ret);
		return ret;
	}
	bool atomicity() const{
		int ret;
		MPI_File_get_atomicity(impl_, &ret);
		return ret;
	}
	void atomicity(bool flag){MPI_File_set_atomicity(impl_, flag);}
	MPI_Offset byte_offset(MPI_Offset offset){
		MPI_Offset disp;
		MPI_File_get_byte_offset(impl_, offset, &disp);
		return disp;
	}
	info hints() const{
		info ret;
		MPI_File_get_info(impl_, &ret.impl_);
		return ret;
	}
	void hints(info set){MPI_File_set_info(impl_, set.impl_);}
	MPI_Offset position() const{
		MPI_Offset offset;
		MPI_File_get_position(impl_, &offset);
		return offset;
	}
	MPI_Offset position_shared() const{
		MPI_Offset offset;
		MPI_File_get_position_shared(impl_, &offset);
		return offset;
	}
	MPI_Offset size() const{
		MPI_Offset ret;
		MPI_File_get_size(impl_, &ret);
		return ret;
	}
	MPI_Aint extent(boost::mpi3::type const& t) const{
		MPI_Aint ret;
		MPI_File_get_type_extent(impl_, t.impl_, &ret);
		return ret;
	}
	void view(); // int MPI_File_get_view
	template<class ContiguousIterator, class Size, class value_type = typename std::iterator_traits<ContiguousIterator>::value_type, class datatype = detail::datatype<value_type>>
	auto read_n(ContiguousIterator O, Size count) const{
		status s;
		MPI_File_read(impl_, std::addressof(*O), count, datatype{}, &s.impl_);
		return O + count;
	}
	template<class ContiguousIterator, class Size, class value_type = typename std::iterator_traits<ContiguousIterator>::value_type, class datatype = detail::datatype<value_type>>
	auto read_all_n(ContiguousIterator O, Size count) const{
		status s;
		MPI_File_read_all(impl_, std::addressof(*O), count, datatype{}, &s.impl_);
		return O + count;
	}
	// MPI_File_read_all_begin
	// MPI_File_read_all_end
	template<class ContiguousIterator, class Size, class value_type = typename std::iterator_traits<ContiguousIterator>::value_type, class datatype = detail::datatype<value_type>>
	auto read_at_n(MPI_Offset offset, ContiguousIterator O, Size count) const{
		status s;
		MPI_File_read_at(offset, impl_, std::addressof(*O), count, datatype{}, &s.impl_);
		return O + count;
	}
	template<class ContiguousIterator, class Size, class value_type = typename std::iterator_traits<ContiguousIterator>::value_type, class datatype = detail::datatype<value_type>>
	auto read_at_all_n(MPI_Offset offset, ContiguousIterator O, Size count) const{
		status s;
		MPI_File_read_at_all(offset, impl_, std::addressof(*O), count, datatype{}, &s.impl_);
		return O + count;
	}
	// MPI_File_read_at_all_begin
	// MPI_File_read_at_all_end
	template<class ContiguousIterator, class Size, class value_type = typename std::iterator_traits<ContiguousIterator>::value_type, class datatype = detail::datatype<value_type>>
	auto read_ordered_n(ContiguousIterator O, Size count) const{
		status s;
		MPI_File_read_ordered(impl_, std::addressof(*O), count, datatype{}, &s.impl_);
		return O + count;
	}
	// MPI_File_read_ordered_begin
	// MPI_File_read_ordered_end
	template<class ContiguousIterator, class Size, class value_type = typename std::iterator_traits<ContiguousIterator>::value_type, class datatype = detail::datatype<value_type>>
	auto read_shared_n(ContiguousIterator O, Size count) const{
		status s;
		MPI_File_read_ordered(impl_, std::addressof(*O), count, datatype{}, &s.impl_);
		return O + count;
	}
	void seek(MPI_Offset offset, int whence){MPI_File_seek(impl_, offset, whence);}
	void seek_set(MPI_Offset offset){seek(offset, MPI_SEEK_SET);}
	void seek_current(MPI_Offset offset){seek(offset, MPI_SEEK_CUR);}
	void seek_end(MPI_Offset offset = 0){seek(offset, MPI_SEEK_END);}

	void seek_shared(MPI_Offset offset, int whence){MPI_File_seek_shared(impl_, offset, whence);}
	//MPI_File_set_errhandler
	template<class ContiguousIterator, class Size, class value_type = typename std::iterator_traits<ContiguousIterator>::value_type, class datatype = detail::datatype<value_type>>
	request iread_n(ContiguousIterator I, Size count);
	template<class ContiguousIterator, class Size, class value_type = typename std::iterator_traits<ContiguousIterator>::value_type, class datatype = detail::datatype<value_type>>
	request iread_at_n(MPI_Offset offset, ContiguousIterator I, Size count);
	template<class ContiguousIterator, class Size, class value_type = typename std::iterator_traits<ContiguousIterator>::value_type, class datatype = detail::datatype<value_type>>
	request iread_shared_n(ContiguousIterator I, Size count);
	template<class ContiguousIterator, class Size, class value_type = typename std::iterator_traits<ContiguousIterator>::value_type, class datatype = detail::datatype<value_type>>
	status write_n(ContiguousIterator I, Size count) const{
		status ret;
		MPI_File_write(impl_, std::addressof(*I), count, datatype{}, &ret.impl_);
		return ret;
	}
	template<class ContiguousIterator, class Size, class value_type = typename std::iterator_traits<ContiguousIterator>::value_type, class datatype = detail::datatype<value_type>>
	status write_all_n(ContiguousIterator I, Size count) const{
		status ret;
		MPI_File_write_all(impl_, std::addressof(*I), count, datatype{}, &ret.impl_);
		return ret;
	}
	template<class ContiguousIterator, class Size, class value_type = typename std::iterator_traits<ContiguousIterator>::value_type, class datatype = detail::datatype<value_type>>
	request iwrite_n(ContiguousIterator I, Size count);
	template<class ContiguousIterator, class Size, class value_type = typename std::iterator_traits<ContiguousIterator>::value_type, class datatype = detail::datatype<value_type>>
	request iwrite_at_n(MPI_Offset offset, ContiguousIterator I, Size count);
	template<class ContiguousIterator, class Size, class value_type = typename std::iterator_traits<ContiguousIterator>::value_type, class datatype = detail::datatype<value_type>>
	request iwrite_shared_n(ContiguousIterator I, Size count);
	void preallocate(MPI_Offset size){MPI_File_preallocate(impl_, size);}
	void resize(MPI_Offset size){MPI_File_set_size(impl_, size);}
	void sync(){MPI_File_sync(impl_);}
};

int ftell(FILE* stream){return stream->position();}
int fseek(FILE* stream, MPI_Offset offset, int whence/*origin*/){
	return MPI_File_seek(stream->impl_, offset, whence);
}

FILE* communicator::fopen(const char* filename, int amode = MPI_MODE_RDWR | MPI_MODE_CREATE){
	FILE* ret;
	MPI_File_open(impl_, filename, amode, MPI_INFO_NULL, &(ret->impl_));
	return ret;
}


}}
#endif

