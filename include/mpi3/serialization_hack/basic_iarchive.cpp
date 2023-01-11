/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// basic_archive.cpp:

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.
// NOLINTBEGIN(altera-id-dependent-backward-branch,altera-unroll-loops,hicpp-use-auto,hicpp-use-equals-default,modernize-use-auto,modernize-use-equals-default,misc-const-correctness)  external code
#include <boost/config.hpp> // msvc 6.0 needs this to suppress warnings

#include <boost/assert.hpp>

#include <cstddef> // size_t, NULL
#include <list>
#include <set>
#include <vector>

#include <boost/config.hpp>  // NOLINT(readability-duplicate-include) third party code
#if defined(BOOST_NO_STDC_NAMESPACE)
namespace std{ 
    using ::size_t; 
} // namespace std
#endif

#include <boost/integer_traits.hpp>

#define BOOST_ARCHIVE_SOURCE
// include this to prevent linker errors when the
// same modules are marked export and import.
#define BOOST_SERIALIZATION_SOURCE
#include <boost/serialization/config.hpp>

#include <boost/serialization/state_saver.hpp>
#include <boost/serialization/throw_exception.hpp>
#include <boost/serialization/tracking.hpp>

#include <boost/archive/archive_exception.hpp>
#include <boost/archive/basic_archive.hpp>
#include <boost/archive/detail/basic_iarchive.hpp>
#include <boost/archive/detail/basic_iserializer.hpp>
#include <boost/archive/detail/basic_pointer_iserializer.hpp>
#include <boost/archive/detail/decl.hpp>

#include <boost/archive/detail/auto_link_archive.hpp>

using namespace boost::serialization;  // NOLINT(google-build-using-namespace,google-global-names-in-headers) external code TODO(correaa) avoid global namespace

namespace boost {
namespace archive {
namespace detail {

class basic_iarchive_impl {
    friend class basic_iarchive;
    library_version_type m_archive_library_version;
    unsigned int m_flags;

    //////////////////////////////////////////////////////////////////////
    // information about each serialized object loaded
    // indexed on object_id
    struct aobject
    {
        void * address;  // NOLINT(modernize-use-default-member-init,misc-non-private-member-variables-in-classes) external code
        bool loaded_as_pointer;  // NOLINT(modernize-use-default-member-init,misc-non-private-member-variables-in-classes) external code
        class_id_type class_id;  // NOLINT(misc-non-private-member-variables-in-classes) external code
        aobject(
            void *a,
            class_id_type class_id_  // NOLINT(performance-unnecessary-value-param) external code
        ) :
            address(a),
            loaded_as_pointer(false),
            class_id(class_id_)
        {}
        aobject() : 
            address(nullptr),
            loaded_as_pointer(false),
            class_id(-2) 
        {}
    };
    typedef std::vector<aobject> object_id_vector_type;  // NOLINT(modernize-use-using) external code
    object_id_vector_type object_id_vector;

    //////////////////////////////////////////////////////////////////////
    // used to implement the reset_object_address operation.
    struct moveable_objects {
        object_id_type start;  // NOLINT(misc-non-private-member-variables-in-classes) external code
        object_id_type end;  // NOLINT(misc-non-private-member-variables-in-classes) external code
        object_id_type recent;  // NOLINT(misc-non-private-member-variables-in-classes) external code
        bool is_pointer;  // NOLINT(modernize-use-default-member-init,misc-non-private-member-variables-in-classes) external code
        moveable_objects() :
            start(0),
            end(0),
            recent(0),
            is_pointer(false)
        {}
    } m_moveable_objects;

    void reset_object_address(
        const void * new_address, 
        const void *old_address
    );

    //////////////////////////////////////////////////////////////////////
    // used by load object to look up class id given basic_serializer
    struct cobject_type  // NOLINT(cppcoreguidelines-special-member-functions,hicpp-special-member-functions) external code
    {
        const basic_iserializer * m_bis;  // NOLINT(misc-non-private-member-variables-in-classes) external code
        const class_id_type m_class_id;  // NOLINT(misc-non-private-member-variables-in-classes) external code
        cobject_type(
            std::size_t class_id,
            const basic_iserializer & bis
        ) : 
            m_bis(& bis),
            m_class_id(class_id)
        {}
        cobject_type(const cobject_type & rhs) :  // NOLINT(hicpp-use-equals-default,modernize-use-equals-default) external code
            m_bis(rhs.m_bis),
            m_class_id(rhs.m_class_id)
        {}
        // the following cannot be defined because of the const
        // member.  This will generate a link error if an attempt
        // is made to assign.  This should never be necessary
        cobject_type & operator=(const cobject_type & rhs);
        bool operator<(const cobject_type &rhs) const
        {
            return *m_bis < *(rhs.m_bis);
        }
    };
    typedef std::set<cobject_type> cobject_info_set_type;  // NOLINT(modernize-use-using) external code
    cobject_info_set_type cobject_info_set;

    //////////////////////////////////////////////////////////////////////
    // information about each serialized class indexed on class_id
    class cobject_id  // NOLINT(cppcoreguidelines-special-member-functions,hicpp-special-member-functions) external code
    {
    public:
        cobject_id & operator=(const cobject_id & rhs){  // NOLINT(hicpp-use-equals-default,modernize-use-equals-default,bugprone-unhandled-self-assignment,cert-oop54-cpp) external code
            bis_ptr = rhs.bis_ptr;
            bpis_ptr = rhs.bpis_ptr;
            file_version = rhs.file_version;
            tracking_level = rhs.tracking_level;
            initialized = rhs.initialized;
            return *this;
        }
        const basic_iserializer * bis_ptr;           // NOLINT(misc-non-private-member-variables-in-classes) third party code
        const basic_pointer_iserializer * bpis_ptr;  // NOLINT(misc-non-private-member-variables-in-classes,modernize-use-default-member-init) third party code
        version_type file_version;                   // NOLINT(misc-non-private-member-variables-in-classes) third party code
        tracking_type tracking_level;                // NOLINT(misc-non-private-member-variables-in-classes) third party code
        bool initialized;                            // NOLINT(misc-non-private-member-variables-in-classes,modernize-use-default-member-init) third party code

        explicit cobject_id(const basic_iserializer & bis_) :
            bis_ptr(& bis_),
            bpis_ptr(nullptr),
            file_version(0),
            tracking_level(track_never),  // NOLINT(readability-implicit-bool-conversion) external code
            initialized(false)
        {}
        cobject_id(const cobject_id &rhs):  // NOLINT(hicpp-use-equals-default,modernize-use-equals-default) external code
            bis_ptr(rhs.bis_ptr),
            bpis_ptr(rhs.bpis_ptr),
            file_version(rhs.file_version),
            tracking_level(rhs.tracking_level),
            initialized(rhs.initialized)
        {}
    };
    typedef std::vector<cobject_id> cobject_id_vector_type;  // NOLINT(modernize-use-using) external code
    cobject_id_vector_type cobject_id_vector;

    //////////////////////////////////////////////////////////////////////
    // address of the most recent object serialized as a poiner
    // whose data itself is now pending serialization
    struct pending {
        void * object;  // NOLINT(misc-non-private-member-variables-in-classes,modernize-use-default-member-init) external code
        const basic_iserializer * bis;  // NOLINT(misc-non-private-member-variables-in-classes,modernize-use-default-member-init) external code
        version_type version;  // NOLINT(misc-non-private-member-variables-in-classes) external code
        pending() :
            object(nullptr),
            bis(nullptr),
            version(0)
        {}
    } m_pending;

    explicit basic_iarchive_impl(unsigned int flags) :
        m_archive_library_version(BOOST_ARCHIVE_VERSION()),
        m_flags(flags)
    {}
    void set_library_version(library_version_type archive_library_version){  // NOLINT(performance-unnecessary-value-param) external code
        m_archive_library_version = archive_library_version;
    }
    bool
    track(
        basic_iarchive & ar,
        void * & t
    );
    void
    load_preamble(
        basic_iarchive & ar,
        cobject_id & co
    );
    class_id_type register_type(
        const basic_iserializer & bis
    );

    // redirect through virtual functions to load functions for this archive
    template<class T>
    void load(basic_iarchive & ar, T & t){
        ar.vload(t);
    }

//public:
    void
    next_object_pointer(void * t){
        m_pending.object = t;
    }
    void delete_created_pointers();
    class_id_type register_type(
        const basic_pointer_iserializer & bpis
    );
    void load_object(
        basic_iarchive & ar,
        void * t,
        const basic_iserializer & bis
    );
    const basic_pointer_iserializer * load_pointer(
        basic_iarchive & ar,
        void * & t, 
        const basic_pointer_iserializer * bpis,
        const basic_pointer_iserializer * (*finder)(
            const boost::serialization::extended_type_info & type
        )
    );
};

inline void 
basic_iarchive_impl::reset_object_address(
    void const * const new_address,  // NOLINT(bugprone-easily-swappable-parameters) external code
    void const * const old_address
){
    if(m_moveable_objects.is_pointer) {
        return; }

    // this code handles a couple of situations.
    // a) where reset_object_address is applied to an untracked object.
    //    In such a case the call is really superfluous and its really an
    //    an error.  But we don't have access to the types here so we can't
    //    know that.  However, this code will effectively turn this situation
    //    into a no-op and every thing will work fine - albeat with a small
    //    execution time penalty.
    // b) where the call to reset_object_address doesn't immediatly follow
    //    the << operator to which it corresponds.  This would be a bad idea
    //    but the code may work anyway.  Naturally, a bad practice on the part
    //    of the programmer but we can't detect it - as above.  So maybe we
    //    can save a few more people from themselves as above.
    object_id_type i = m_moveable_objects.recent;
    for(; i < m_moveable_objects.end; ++i){
        if(old_address == object_id_vector[i].address) {
            break; }
    }
    for(; i < m_moveable_objects.end; ++i){
        void const * const this_address = object_id_vector[i].address;
        // calculate displacement from this level
        // warning - pointer arithmetic on void * is in herently non-portable
        // but expected to work on all platforms in current usage
        if(this_address > old_address){
            std::size_t member_displacement
                = reinterpret_cast<std::size_t>(this_address)  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) external
                - reinterpret_cast<std::size_t>(old_address);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) external
            object_id_vector[i].address = reinterpret_cast<void *>(  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,performance-no-int-to-ptr) external
                reinterpret_cast<std::size_t>(new_address) + member_displacement  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) external
            );
        }
        else{
            std::size_t member_displacement
                = reinterpret_cast<std::size_t>(old_address)  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) external
                - reinterpret_cast<std::size_t>(this_address);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) external
            object_id_vector[i].address = reinterpret_cast<void *>(  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,performance-no-int-to-ptr) external
                reinterpret_cast<std::size_t>(new_address) - member_displacement  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) external code
            );
       }
    }
}

inline void 
basic_iarchive_impl::delete_created_pointers()
{
    object_id_vector_type::iterator i;
    for(
        i = object_id_vector.begin();
        i != object_id_vector.end(); 
        ++i
    ){
        if(i->loaded_as_pointer){
            // borland complains without this minor hack
            const int j = i->class_id;
            const cobject_id & co = cobject_id_vector[static_cast<std::size_t>(j)];
            //const cobject_id & co = cobject_id_vector[i->class_id];
            // with the appropriate input serializer, 
            // delete the indicated object
            co.bis_ptr->destroy(i->address);
        }
    }
}

inline class_id_type
basic_iarchive_impl::register_type(
    const basic_iserializer & bis
){
    class_id_type cid(cobject_info_set.size());
    cobject_type co(static_cast<std::size_t>(cid), bis);
    std::pair<cobject_info_set_type::const_iterator, bool>
        result = cobject_info_set.insert(co);

    if(result.second){
        cobject_id_vector.push_back(cobject_id(bis));
        BOOST_ASSERT(cobject_info_set.size() == cobject_id_vector.size());
    }
    cid = result.first->m_class_id;
    // borland complains without this minor hack
    const int tid = cid;
    cobject_id & coid = cobject_id_vector[static_cast<std::size_t>(tid)];
    coid.bpis_ptr = bis.get_bpis_ptr();
    return cid;
}

void
basic_iarchive_impl::load_preamble(
    basic_iarchive & ar,
    cobject_id & co
){
    if(! co.initialized){
        if(co.bis_ptr->class_info()){
            class_id_optional_type cid(class_id_type(0));
            load(ar, cid);    // to be thrown away
            load(ar, co.tracking_level);
            load(ar, co.file_version);
        }
        else{
            // override tracking with indicator from class information
            co.tracking_level = co.bis_ptr->tracking(m_flags);
            co.file_version = version_type(  // NOLINT(google-readability-casting) third pary code
                co.bis_ptr->version()
            );
        }
        co.initialized = true;
    }
}

bool
basic_iarchive_impl::track(
    basic_iarchive & ar,
    void * & t
){
    object_id_type oid;
    load(ar, oid);

    // if its a reference to a old object
    if(object_id_type(object_id_vector.size()) > oid){
        // we're done
        t = object_id_vector[oid].address;
        return false;
    }
    return true;
}

inline void
basic_iarchive_impl::load_object(
    basic_iarchive & ar,
    void * t,
    const basic_iserializer & bis
){
    m_moveable_objects.is_pointer = false;
    serialization::state_saver<bool> ss_is_pointer(m_moveable_objects.is_pointer);
    // if its been serialized through a pointer and the preamble's been done
    if(t == m_pending.object && & bis == m_pending.bis){
        // read data
        (bis.load_object_data)(ar, t, m_pending.version);
        return;
    }

    const class_id_type cid = register_type(bis);
    const int i = cid;
    cobject_id & co = cobject_id_vector[static_cast<std::size_t>(i)];

    load_preamble(ar, co);

    // save the current move stack position in case we want to truncate it
    boost::serialization::state_saver<object_id_type> ss_start(m_moveable_objects.start);

    // note: extra line used to evade borland issue
    const bool tracking = co.tracking_level;

    object_id_type this_id;
    m_moveable_objects.start =
    this_id = object_id_type(object_id_vector.size());

    // if we tracked this object when the archive was saved
    if(tracking){ 
        // if it was already read
        if(!track(ar, t)) {
            // we're done
            return; }
        // add a new enty into the tracking list
        object_id_vector.push_back(aobject(t, cid));
        // and add an entry for this object
        m_moveable_objects.end = object_id_type(object_id_vector.size());
    }
    // read data
    (bis.load_object_data)(ar, t, co.file_version);
    m_moveable_objects.recent = this_id;
}

inline const basic_pointer_iserializer *
basic_iarchive_impl::load_pointer(
    basic_iarchive &ar,
    void * & t,
    const basic_pointer_iserializer * bpis_ptr,
    const basic_pointer_iserializer * (*finder)(
        const boost::serialization::extended_type_info & type_
    )
){
    m_moveable_objects.is_pointer = true;
    serialization::state_saver<bool> w(m_moveable_objects.is_pointer);

    class_id_type cid;
    load(ar, cid);

#if BOOST_VERSION < 107400
    if(NULL_POINTER_TAG == cid){
        t = NULL;
        return bpis_ptr;
    }
#else // this case is taken from https://github.com/boostorg/serialization/blob/develop/src/basic_iarchive.cpp#L430-L433
    if(BOOST_SERIALIZATION_NULL_POINTER_TAG == cid){
        t = nullptr;
        return bpis_ptr;
    }
#endif

    // if its a new class type - i.e. never been registered
    if(class_id_type(cobject_info_set.size()) <= cid){
        // if its either abstract
        if(nullptr == bpis_ptr
        // or polymorphic
        || bpis_ptr->get_basic_serializer().is_polymorphic()){
            // is must have been exported
            char key[BOOST_SERIALIZATION_MAX_KEY_SIZE];  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) external code
            class_name_type class_name(key);
            load(ar, class_name);
            // if it has a class name
            const serialization::extended_type_info *eti = nullptr;
            if(0 != key[0]) {
                eti = serialization::extended_type_info::find(key); }
            if(nullptr == eti) {
                boost::serialization::throw_exception(
                    archive_exception(archive_exception::unregistered_class)
                ); }
            bpis_ptr = (*finder)(*eti);  // cppcheck-suppress [nullPointerRedundantCheck,nullPointer] ; legacy code
        }
        BOOST_ASSERT(nullptr != bpis_ptr);
        // class_id_type new_cid = register_type(bpis_ptr->get_basic_serializer());
        BOOST_VERIFY(register_type(bpis_ptr->get_basic_serializer()) == cid);
        int i = cid;
        cobject_id_vector[static_cast<std::size_t>(i)].bpis_ptr = bpis_ptr;
    }
    int i = cid;
    cobject_id & co = cobject_id_vector[static_cast<std::size_t>(i)];
    bpis_ptr = co.bpis_ptr;

    load_preamble(ar, co);

    // extra line to evade borland issue
    const bool tracking = co.tracking_level;
    // if we're tracking and the pointer has already been read
    if(tracking && ! track(ar, t)) {
        // we're done
        return bpis_ptr; }

    // save state
    serialization::state_saver<object_id_type> w_start(m_moveable_objects.start);

    // allocate space on the heap for the object - to be constructed later
    t = bpis_ptr->heap_allocation();
    BOOST_ASSERT(nullptr != t);

    if(! tracking){
        bpis_ptr->load_object_ptr(ar, t, co.file_version);
    }
    else{
        serialization::state_saver<void *> x(m_pending.object);
        serialization::state_saver<const basic_iserializer *> y(m_pending.bis);
        serialization::state_saver<version_type> z(m_pending.version);

        m_pending.bis = & bpis_ptr->get_basic_serializer();
        m_pending.version = co.file_version;

        // predict next object id to be created
        const unsigned int ui = static_cast<unsigned int>(object_id_vector.size());

        serialization::state_saver<object_id_type> w_end(m_moveable_objects.end);

        
        // add to list of serialized objects so that we can properly handle
        // cyclic strucures
        object_id_vector.push_back(aobject(t, cid));

        // remember that that the address of these elements could change
        // when we make another call so don't use the address
        bpis_ptr->load_object_ptr(
            ar,
            t,
            m_pending.version
        );
        object_id_vector[ui].loaded_as_pointer = true;
    }

    return bpis_ptr;
}

} // namespace detail
} // namespace archive
} // namespace boost

//////////////////////////////////////////////////////////////////////
// implementation of basic_iarchive functions
namespace boost {
namespace archive {
namespace detail {

BOOST_ARCHIVE_DECL void
basic_iarchive::next_object_pointer(void *t){
    pimpl->next_object_pointer(t);
}

BOOST_ARCHIVE_DECL
basic_iarchive::basic_iarchive(unsigned int flags) : 
    pimpl(new basic_iarchive_impl(flags))
{}

BOOST_ARCHIVE_DECL
basic_iarchive::~basic_iarchive()
{}

BOOST_ARCHIVE_DECL void
basic_iarchive::set_library_version(library_version_type archive_library_version){  // NOLINT(performance-unnecessary-value-param)
    pimpl->set_library_version(archive_library_version);
}

BOOST_ARCHIVE_DECL void
basic_iarchive::reset_object_address(
    const void * new_address, 
    const void * old_address
){
    pimpl->reset_object_address(new_address, old_address);
}

BOOST_ARCHIVE_DECL void
basic_iarchive::load_object(
    void *t, 
    const basic_iserializer & bis
){
    pimpl->load_object(*this, t, bis);
}

// load a pointer object
BOOST_ARCHIVE_DECL const basic_pointer_iserializer *
basic_iarchive::load_pointer(
    void * &t, 
    const basic_pointer_iserializer * bpis_ptr,
    const basic_pointer_iserializer * (*finder)(
        const boost::serialization::extended_type_info & type_
    )

){
    return pimpl->load_pointer(*this, t, bpis_ptr, finder);
}

BOOST_ARCHIVE_DECL void
basic_iarchive::register_basic_serializer(const basic_iserializer & bis){
    pimpl->register_type(bis);
}

BOOST_ARCHIVE_DECL void
basic_iarchive::delete_created_pointers()
{
    pimpl->delete_created_pointers();
}

BOOST_ARCHIVE_DECL boost::archive::library_version_type
basic_iarchive::get_library_version() const{
    return pimpl->m_archive_library_version;
}

BOOST_ARCHIVE_DECL unsigned int
basic_iarchive::get_flags() const{
    return pimpl->m_flags;
}

} // namespace detail
} // namespace archive
} // namespace boost
// NOLINTEND(altera-id-dependent-backward-branch,altera-unroll-loops,hicpp-use-auto,hicpp-use-equals-default,modernize-use-auto,modernize-use-equals-default,misc-const-correctness)  external code
