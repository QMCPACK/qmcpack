
#ifndef BP_OPEN_H
#define BP_OPEN_H

#include "JabberWocky.h"
#include "utils.h"

#define MINIFOOTER_SIZE 28
#define ADIOS_VERSION_NUM_MASK                       0x000000FF
#define ADIOS_VERSION_BP_FORMAT                      2
#define BYTE_ALIGN 8

void swap_ptr(void * data, int size);

void swap_16_ptr(void *data);

void swap_32_ptr(void *data);

void swap_64_ptr(void *data);

void swap_128_ptr(void *data);

#define swap_16(data) swap_16_ptr(&(data))

#define swap_32(data) swap_32_ptr(&(data))

#define swap_64(data) swap_64_ptr(&(data))

#define swap_128(data) swap_128_ptr(&(data))

#define BUFREAD8(b,var)  var = (uint8_t) *(b->buff + b->offset); \
                                                        b->offset += 1;

#define BUFREAD16(b,var) var = *(uint16_t *) (b->buff + b->offset); \
                                                        if (b->change_endianness == adios_flag_yes) \
                             swap_16(var); \
                         b->offset += 2;

#define BUFREAD32(b,var) var = *(uint32_t *) (b->buff + b->offset); \
                                                        if (b->change_endianness == adios_flag_yes) \
                             swap_32(var); \
                         b->offset += 4;

#define BUFREAD64(b,var) var = *(uint64_t *) (b->buff + b->offset); \
                                                        if (b->change_endianness == adios_flag_yes) \
                             swap_64(var); \
                         b->offset += 8;


struct BP_file_handle
{
    uint32_t file_index;
    MPI_File fh;
    struct BP_file_handle * next;
};

typedef struct BP_file_handle BP_file_handle_list;

struct bp_minifooter {
    uint64_t time_steps;  /* = fh->tidx_stop - fh->tidx_start + 1 */
    uint64_t pgs_count;
    uint64_t pgs_length;
    uint32_t vars_count;
    uint32_t attrs_count;
    uint64_t vars_length;
    uint64_t attrs_length;
    uint64_t pgs_index_offset;
    uint64_t vars_index_offset;
    uint64_t attrs_index_offset;
    uint32_t version;
    uint32_t change_endianness; // = enum ADIOS_FLAG, 0: unknown!, adios_flag_yes or adios_flag_no
    uint64_t file_size;
} __attribute__((__packed__));

struct bp_index_pg_struct_v1
{
    char * group_name;
    enum ADIOS_FLAG adios_host_language_fortran;
    uint32_t process_id;
    char * time_index_name;
    uint32_t time_index;
    uint64_t offset_in_file;

    struct bp_index_pg_struct_v1 * next;
};

struct adios_index_var_struct_v1
{
    uint32_t id;
    char * group_name;
    char * var_name;
    char * var_path;
    enum ADIOS_DATATYPES type;

    uint64_t characteristics_count;
    uint64_t characteristics_allocated;

    struct adios_index_characteristic_struct_v1 * characteristics;

    struct adios_index_var_struct_v1 * next;
};

struct adios_index_attribute_struct_v1
{
    uint32_t id;
    char * group_name;
    char * attr_name;
    char * attr_path;
    enum ADIOS_DATATYPES type;

    uint64_t characteristics_count;
    uint64_t characteristics_allocated;

    struct adios_index_characteristic_struct_v1 * characteristics;

    struct adios_index_attribute_struct_v1 * next;
};

struct BP_GROUP_VAR {
    uint16_t group_count;
    uint16_t group_id;
    char ** namelist;
    uint32_t *** time_index;
    uint64_t * pg_offsets;
    char ** var_namelist;
    uint32_t * var_counts_per_group;
    uint64_t ** var_offsets;
};

struct BP_GROUP_ATTR {
    uint16_t group_count;
    uint16_t group_id;
    char ** namelist;
    char ** attr_namelist;
    uint32_t * attr_counts_per_group;
    uint64_t ** attr_offsets;
};

struct BP_GROUP {
    uint16_t group_id;
    uint16_t vars_offset; // start of variables belonging to this group in the list of variables from all groups; old read API used this
    uint32_t vars_count;
    uint16_t attrs_offset; // old read API, this group's attributes in the list of all groups' attrs 
    uint32_t attrs_count;
    struct BP_FILE * fh;
    struct adios_index_var_struct_v1 * vars_root;  /* pointer into the list of BP_FILE.vars_root */
    struct adios_index_attribute_struct_v1 * attrs_root;
};

typedef struct BP_FILE {
    MPI_File mpi_fh;
    char * fname; // Main file name is needed to calculate subfile names
    BP_file_handle_list * sfh; // This list links all the subfiles handle together
    MPI_Comm comm;
    struct adios_bp_buffer_struct_v1 * b;
    struct bp_index_pg_struct_v1 * pgs_root;
    struct adios_index_var_struct_v1 * vars_root;
    struct adios_index_attribute_struct_v1 * attrs_root;
    struct adios_index_var_struct_v1 ** vars_table; // To speed up vars_root lookup. Q. Liu, 12-2013.
    struct bp_minifooter mfooter;
    struct BP_GROUP_VAR * gvar_h;
    struct BP_GROUP_ATTR * gattr_h;
    uint32_t tidx_start;
    uint32_t tidx_stop;
    void * priv;
} BP_FILE;

struct adios_bp_buffer_struct_v1
{
    int f;             // the file handle
    uint64_t file_size;
    uint32_t version;

    char * allocated_buff_ptr;  // initial alloc for aligning on 8-byte boundary

    char * buff;
    uint64_t length;
    uint64_t offset;   // buffer_offset

    enum ADIOS_FLAG change_endianness;

    off_t file_offset;
    uint64_t end_of_pgs;          // where the last process group ends
    uint64_t pg_index_offset;     // process groups index starts
    uint64_t pg_size;             // process groups index size
    uint64_t vars_index_offset;   // vars index start
    uint64_t vars_size;           // vars index size
    uint64_t attrs_index_offset;  // attributes index start
    uint64_t attrs_size;          // attributes index size

    uint64_t read_pg_offset;
    uint64_t read_pg_size;
};


int bp_open(const char *fname, MPI_Comm comm, struct BP_FILE *fh);
void adios_buffer_struct_init (struct adios_bp_buffer_struct_v1 * b);
int bp_read_open (const char * filename, MPI_Comm comm,struct BP_FILE * fh);
int bp_read_minifooter (struct BP_FILE * bp_struct);
void bp_alloc_aligned (struct adios_bp_buffer_struct_v1 * b, uint64_t size);
void bp_realloc_aligned (struct adios_bp_buffer_struct_v1 * b, uint64_t size);
int adios_parse_version (struct adios_bp_buffer_struct_v1 * b, uint32_t * version);
#endif
