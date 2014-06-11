
#include "include/bp_open.h"
#include "include/utils.h"

int bp_open (const char * fname,
             MPI_Comm comm,
             struct BP_FILE * fh)
{
    int rank;
    uint64_t header_size;

    MPI_Comm_rank (comm, &rank);

    adios_buffer_struct_init (fh->b);

    if (bp_read_open (fname, comm, fh))
    {
        return -1;
    }

    /* Only rank 0 reads the footer and it broadcasts to all other processes */
    if (rank == 0)
    {
        if (bp_read_minifooter (fh))
        {
            return -1;
        }
    }

    /* Broadcast to all other processors */
    MPI_Bcast (&fh->mfooter, sizeof (struct bp_minifooter), MPI_BYTE, 0, comm);

    header_size = fh->mfooter.file_size-fh->mfooter.pgs_index_offset;

    if (rank != 0)
    {
        if (!fh->b->buff)
        {
            bp_alloc_aligned (fh->b, header_size);
            assert (fh->b->buff);

            memset (fh->b->buff, 0, header_size);
            fh->b->offset = 0;
        }
    }

    MPI_Barrier (comm);
    MPI_Bcast (fh->b->buff, fh->mfooter.file_size-fh->mfooter.pgs_index_offset, MPI_BYTE, 0, comm);

    /* Everyone parses the index on its own */
    //bp_parse_pgs (fh);
    //bp_parse_vars (fh);
    //bp_parse_attrs (fh);

    return 0;
}

void adios_buffer_struct_init (struct adios_bp_buffer_struct_v1 * b)
{
    b->f = -1;
    b->allocated_buff_ptr = 0;
    b->buff = 0;
    b->length = 0;
    b->change_endianness = adios_flag_unknown;
    b->version = 0;
    b->offset = 0;
    b->end_of_pgs = 0;
    b->pg_index_offset = 0;
    b->pg_size = 0;
    b->vars_index_offset = 0;
    b->vars_size = 0;
    b->file_size = 0;
    b->read_pg_offset = 0;
    b->read_pg_size = 0;
}

int bp_read_open (const char * filename,
          MPI_Comm comm,
          struct BP_FILE * fh)
{
    int  err;
    int  rank;

    MPI_Comm_rank (comm, &rank);

    // variable definition
    MPI_Offset  file_size;

    // open a file by the multiple processors within the same
    // communicator
    err = MPI_File_open (comm, (char *) filename, MPI_MODE_RDONLY,
            (MPI_Info) MPI_INFO_NULL, &(fh->mpi_fh));
    if (err != MPI_SUCCESS) {
        char e [MPI_MAX_ERROR_STRING];
        int len = 0;
        memset (e, 0, MPI_MAX_ERROR_STRING);
        MPI_Error_string (err, e, &len);
        adios_error (err_file_open_error, "MPI open failed for %s: '%s'\n", filename, e);
        return adios_flag_no;
    }

    MPI_File_get_size (fh->mpi_fh, &file_size);
    fh->b->file_size = file_size;
    fh->mfooter.file_size = file_size;

    return 0;
}

int bp_read_minifooter (struct BP_FILE * bp_struct)
{
    struct adios_bp_buffer_struct_v1 * b = bp_struct->b;
    struct bp_minifooter * mh = &bp_struct->mfooter;
    uint64_t attrs_end = b->file_size - MINIFOOTER_SIZE;
    int r;

    MPI_Status status;

    if (!b->buff) {
        bp_alloc_aligned (b, MINIFOOTER_SIZE);
        if (!b->buff) {
            adios_error (err_no_memory, "could not allocate %d bytes\n", MINIFOOTER_SIZE);
            return 1;
        }
        memset (b->buff, 0, MINIFOOTER_SIZE);
        b->offset = 0;
    }
    MPI_File_seek (bp_struct->mpi_fh, (MPI_Offset) attrs_end, MPI_SEEK_SET);
    MPI_File_read (bp_struct->mpi_fh, b->buff, MINIFOOTER_SIZE, MPI_BYTE, &status);

    /*memset (&mh->pgs_index_offset, 0, MINIFOOTER_SIZE);
    memcpy (&mh->pgs_index_offset, b->buff, MINIFOOTER_SIZE);*/

    /* get version id. Needs the bp->offset be pointing to the last 4 bytes of the buffer,
       It also sets b->change_endianness */
    /* Note that b is not sent over to processes, only the minifooter and then b->buff (the footer) */
    b->offset = MINIFOOTER_SIZE - 4;
    adios_parse_version (b, &mh->version);
    mh->change_endianness = b->change_endianness;

    // validity check
    if ((mh->version & ADIOS_VERSION_NUM_MASK) > ADIOS_VERSION_BP_FORMAT) {
        adios_error (err_file_open_error,
           "Invalid BP file detected. Format version of file seems to be %d, "
           "which is greater than the highest supported version %d.\n",
           (mh->version & ADIOS_VERSION_NUM_MASK), ADIOS_VERSION_BP_FORMAT);
        return 1;
    }

    b->offset = 0; // reset offset to beginning

    BUFREAD64(b, b->pg_index_offset)
    mh->pgs_index_offset = b->pg_index_offset;
    // validity check  
    if (b->pg_index_offset > b->file_size) {
        adios_error (err_file_open_error,
                "Invalid BP file detected. PG index offset (%lld) > file size (%lld)\n",
                b->pg_index_offset, b->file_size);
        return 1;
    }

    BUFREAD64(b, b->vars_index_offset)
    mh->vars_index_offset = b->vars_index_offset;
    // validity check  
    if (b->vars_index_offset > b->file_size) {
        adios_error (err_file_open_error,
                "Invalid BP file detected. Variable index offset (%lld) > file size (%lld)\n",
                b->vars_index_offset, b->file_size);
        return 1;
    }
    if (b->vars_index_offset < b->pg_index_offset) {
        adios_error (err_file_open_error,
                "Invalid BP file detected. Variable index offset (%lld) < PG index offset (%lld)\n",
                b->vars_index_offset, b->pg_index_offset);
        return 1;
    }


    BUFREAD64(b, b->attrs_index_offset)
    mh->attrs_index_offset = b->attrs_index_offset;
    // validity check  
    if (b->attrs_index_offset > b->file_size) {
        adios_error (err_file_open_error,
                "Invalid BP file detected. Attribute index offset (%lld) > file size (%lld)\n",
                b->attrs_index_offset, b->file_size);
        return 1;
    }
    if (b->attrs_index_offset < b->vars_index_offset) {
        adios_error (err_file_open_error,
                "Invalid BP file detected. Attribute index offset (%lld) < Variable index offset (%lld)\n",
                b->attrs_index_offset, b->vars_index_offset);
        return 1;
    }

    b->end_of_pgs = b->pg_index_offset;
    b->pg_size = b->vars_index_offset - b->pg_index_offset;
    b->vars_size = b->attrs_index_offset - b->vars_index_offset;
    b->attrs_size = attrs_end - b->attrs_index_offset;

    /* Read the whole footer */
    /* FIXME: including the last 28 bytes read already above and it seems that is not processed anymore */
    /* It will be sent to all processes */
    uint64_t footer_size = mh->file_size - mh->pgs_index_offset;
    bp_realloc_aligned (b, footer_size);
    MPI_File_seek (bp_struct->mpi_fh,
                        (MPI_Offset)  mh->pgs_index_offset,
                        MPI_SEEK_SET);
    MPI_File_read (bp_struct->mpi_fh, b->buff, footer_size,
            MPI_BYTE, &status);

    MPI_Get_count (&status, MPI_BYTE, &r);

    // reset the pointer to the beginning of buffer
    b->offset = 0;
    return 0;
}

void bp_alloc_aligned (struct adios_bp_buffer_struct_v1 * b, uint64_t size)
{

    b->allocated_buff_ptr =  (char *)malloc (size + BYTE_ALIGN - 1);
    if (!b->allocated_buff_ptr)
    {
        adios_error ( err_no_memory, "Cannot allocate %llu bytes\n", size);

        b->buff = NULL;
        b->length = 0;

        return;
    }
    uint64_t p = (uint64_t) b->allocated_buff_ptr;
    b->buff = (char *) ((p + BYTE_ALIGN - 1) & ~(BYTE_ALIGN - 1));
    b->length = size;
}

void bp_realloc_aligned (struct adios_bp_buffer_struct_v1 * b
                            ,uint64_t size
                            )
{
    b->allocated_buff_ptr = (char *)realloc (b->allocated_buff_ptr
                                    ,size + BYTE_ALIGN - 1
                                    );
    if (!b->allocated_buff_ptr)
    {
        adios_error ( err_no_memory, "Cannot allocate %llu bytes\n", size);

        b->buff = NULL;
        b->length = 0;

        return;
    }
    uint64_t p = (uint64_t) b->allocated_buff_ptr;
    b->buff = (char *) ((p + BYTE_ALIGN - 1) & ~(BYTE_ALIGN - 1));
    b->length = size;
}


int adios_parse_version (struct adios_bp_buffer_struct_v1 * b,
                        uint32_t * version
                        )
{
    // if high bit set, big endian
    uint32_t test = 1;

    if (b->length < 4)
    {
        adios_error(err_invalid_buffer_version, "adios_parse_version requires"
                "a buffer of at least "
                "4 bytes.  Only %llu were provided\n", b->length);
        return 1;
    }

    *version = ntohl (*(uint32_t *) (b->buff + b->offset));
    char *v = (char *) version;
    if (   (*v && !*(char *) &test)       // both writer and this machine are big endian
            || (!*(v+3) && *(char *) &test)   // both are little endian
       )
    {
        b->change_endianness = adios_flag_no;//no need to change endiannness
    }
    else
    {
        b->change_endianness = adios_flag_yes;
    }

    *version = *version & 0x7fffffff;

    return 0;
}

void swap_16_ptr(void *data)
{
    uint16_t d = *(uint16_t *)data;
    *(uint16_t *)data = d>>8 | d<<8;
}


void swap_32_ptr(void *data)
{
    uint32_t d = *(uint32_t *)data;
    *(uint32_t *)data = ((d&0x000000FF)<<24)
                      + ((d&0x0000FF00)<<8)
                      + ((d&0x00FF0000)>>8)
                      + ((d&0xFF000000)>>24);
}


void swap_64_ptr(void *data)
{
    uint64_t d = *(uint64_t *)data;
    *(uint64_t *)data = ((d&0x00000000000000FF)<<56)
                          + ((d&0x000000000000FF00)<<40)
                          + ((d&0x0000000000FF0000)<<24)
                          + ((d&0x00000000FF000000)<<8)
                          + ((d&0x000000FF00000000LL)>>8)
                          + ((d&0x0000FF0000000000LL)>>24)
                          + ((d&0x00FF000000000000LL)>>40)
                          + ((d&0xFF00000000000000LL)>>56);
}


void swap_128_ptr(void *data)
{
    uint64_t d = *(uint64_t *)data;
    *(uint64_t *)data = ((d&0x00000000000000FF)<<56)
                          + ((d&0x000000000000FF00)<<40)
                          + ((d&0x0000000000FF0000)<<24)
                          + ((d&0x00000000FF000000)<<8)
                          + ((d&0x000000FF00000000LL)>>8)
                          + ((d&0x0000FF0000000000LL)>>24)
                          + ((d&0x00FF000000000000LL)>>40)
                          + ((d&0xFF00000000000000LL)>>56);
    d = *((uint64_t *)data + 1);
    d = ((d&0x00000000000000FF)<<56)
                          + ((d&0x000000000000FF00)<<40)
                          + ((d&0x0000000000FF0000)<<24)
                          + ((d&0x00000000FF000000)<<8)
                          + ((d&0x000000FF00000000LL)>>8)
                          + ((d&0x0000FF0000000000LL)>>24)
                          + ((d&0x00FF000000000000LL)>>40)
                          + ((d&0xFF00000000000000LL)>>56);
    *((uint64_t *)data + 1) = *(uint64_t *)data;
    *(uint64_t *)data = d;
}

