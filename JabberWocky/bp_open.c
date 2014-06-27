
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

    if (bp_read_minifooter (fh))
    {
        return -1;
    }

    //MPI_Bcast (&fh->mfooter, sizeof (struct bp_minifooter), MPI_BYTE, 0, comm);

    header_size = fh->mfooter.file_size-fh->mfooter.pgs_index_offset;
    printf("header %llu\n", header_size);

    //if (rank != 0)
    //{
      //  if (!fh->b->buff)
        //{
          //  bp_alloc_aligned (fh->b, header_size);
            //assert (fh->b->buff);

           // memset (fh->b->buff, 0, header_size);
           // fh->b->offset = 0;
       // }
    //}

    //MPI_Barrier (comm);
    //MPI_Bcast (fh->b->buff, fh->mfooter.file_size-fh->mfooter.pgs_index_offset, MPI_BYTE, 0, comm);

    /* Everyone parses the index on its own */
    bp_parse_pgs (fh);
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
    printf("file size %llu\n", file_size);
    fh->b->file_size = file_size;
    fh->mfooter.file_size = file_size;

    return 0;
}

int bp_read_minifooter (struct BP_FILE * bp_struct)
{
    struct adios_bp_buffer_struct_v1 * b = bp_struct->b;
    struct bp_minifooter * mh = &bp_struct->mfooter;
    uint64_t minifooter_start= b->file_size - MINIFOOTER_SIZE;
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
    MPI_File_seek (bp_struct->mpi_fh, (MPI_Offset) minifooter_start, MPI_SEEK_SET);
    MPI_File_read (bp_struct->mpi_fh, b->buff, MINIFOOTER_SIZE, MPI_BYTE, &status);

    memset (&mh->pgs_index_offset, 0, MINIFOOTER_SIZE);
    memcpy (&mh->pgs_index_offset, b->buff, MINIFOOTER_SIZE);
    printf("pg_index_offset %llu\n", mh->pgs_index_offset);

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
    printf("pg_index_offset %llu\n", mh->pgs_index_offset);
    // validity check  
    if (b->pg_index_offset > b->file_size) {
        adios_error (err_file_open_error,
                "Invalid BP file detected. PG index offset (%lld) > file size (%lld)\n",
                b->pg_index_offset, b->file_size);
        return 1;
    }

    BUFREAD64(b, b->vars_index_offset)
    mh->vars_index_offset = b->vars_index_offset;
    printf("vars_index_offset %llu\n", mh->vars_index_offset);
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
    printf("attrs_index_offset %llu\n", mh->attrs_index_offset);
    // validity check s 
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
    b->attrs_size = minifooter_start - GAP - b->attrs_index_offset;
    printf("end_of_pgs %llu\n", b->end_of_pgs); 
    printf("pg_size %llu\n", b->pg_size); 
    printf("vars_size %llu\n", b->vars_size); 
    printf("attrs_size %llu\n", b->attrs_size); 

    /* Read the whole footer */
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

/****************/
/* Parse GROUPS */
/****************/
int bp_parse_pgs (struct BP_FILE * fh)
{
    struct bp_index_pg_struct_v1 ** root = &(fh->pgs_root); // need the pointer to it to malloc below
    struct adios_bp_buffer_struct_v1 * b = fh->b;
    struct bp_minifooter * mh = &(fh->mfooter);
    uint64_t i;

    /* Note that at this point, many variables of b->* is unset (init'ed to 0).
       It's the minifooter which holds accurate information.
       b holds the footer data from the file */

    b->offset = 0;
    b->change_endianness = (enum ADIOS_FLAG) mh->change_endianness;

    BUFREAD64(b, mh->pgs_count)
    BUFREAD64(b, mh->pgs_length)
    printf(" pgs count %llu pgs_length %llu\n", mh->pgs_count, mh->pgs_length);

    int j;
    uint64_t group_count = 0;
    char ** namelist;
    char fortran_flag;

    namelist = (char **) malloc(sizeof(char *)*mh->pgs_count);
    uint16_t * grpidlist = (uint16_t *) malloc(sizeof(uint16_t)*mh->pgs_count);

    uint32_t tidx_start, tidx_stop; /* Determine first and last timestep in file */

    for (i = 0; i < mh->pgs_count; i++) {
        uint16_t length_of_group;
        namelist[i] = 0;
        // validate remaining length
        BUFREAD16(b, length_of_group)

        if (!*root)
        {
            *root = (struct bp_index_pg_struct_v1 *)
                malloc (sizeof(struct bp_index_pg_struct_v1));
            memset (*root, 0, sizeof(struct bp_index_pg_struct_v1));
            (*root)->next = 0;
        }
        uint16_t length_of_name;

        BUFREAD16(b, length_of_name)
        (*root)->group_name = (char *) malloc (length_of_name + 1);
        (*root)->group_name [length_of_name] = '\0';
        memcpy ((*root)->group_name, b->buff + b->offset, length_of_name);
        b->offset += length_of_name;


        if ( group_count == 0 ) {
            namelist[group_count] = (char *) malloc (length_of_name + 1);
            strcpy (namelist[group_count], (*root)->group_name);
            ++group_count;
            grpidlist[i] = group_count-1;
        }
        else {
            for (j=0; j<group_count; j++) {
                if (!strcmp(namelist[j], (*root)->group_name)) {
                    break;
                }
            }
            if (j==group_count) {
                namelist[group_count] = (char *) malloc (length_of_name + 1);
                strcpy (namelist[group_count], (*root)->group_name);
                ++group_count;
                grpidlist[i] = group_count - 1;
            }
            else
                grpidlist[i] = j;

        }

        BUFREAD8(b, fortran_flag)
        (*root)->adios_host_language_fortran =
            (fortran_flag == 'y' ? adios_flag_yes : adios_flag_no);

        BUFREAD32(b, (*root)->process_id)

        BUFREAD16(b, length_of_name)
        (*root)->time_index_name = (char *) malloc (length_of_name + 1);
        (*root)->time_index_name [length_of_name] = '\0';
        memcpy ((*root)->time_index_name, b->buff + b->offset, length_of_name);
        b->offset += length_of_name;

        BUFREAD32(b, (*root)->time_index)

        BUFREAD64(b, (*root)->offset_in_file)

        if (i == 0)
            tidx_start = (*root)->time_index;
        if (i == mh->pgs_count-1) {
            tidx_stop = (*root)->time_index;
            mh->time_steps = tidx_stop - tidx_start + 1;
        }

        root = &(*root)->next;
    }

    /*
    root = &(fh->pgs_root);
    for (i = 0; i < mh->pgs_count; i++) {
        printf("%d\tpg pid=%d addr=%x next=%x\n",i, (*root)->process_id, *root, (*root)->next);
        root = &(*root)->next;
    }
    */

    uint64_t * pg_offsets = 0;
    uint32_t * pg_pids = 0;
    uint32_t *** time_index = 0;
    pg_offsets = (uint64_t *)
        malloc (sizeof(uint64_t)*mh->pgs_count);
    pg_pids = (uint32_t *)
        malloc (sizeof(uint32_t)*mh->pgs_count);
    // time_index[0]: record which pg to start from per timestep per group
    // time_index[1]: record the # of pgs per timesteps per group
    time_index = (uint32_t ***) malloc (sizeof(uint32_t **)*2);

    for (j=0;j<2;j++) {
        time_index[j] = (uint32_t **)
            malloc (sizeof(uint32_t*)*group_count);
        //printf ("### time_index[%d]=%x  group_count=%d  #pgs=%d #ts=%d\n", j, time_index[j], group_count, mh->pgs_count,  mh->time_steps);
        for (i=0;i<group_count;i++) {
            if (mh->pgs_count < mh->time_steps) {
                /* FIXME: when can this happen?
                   pgs = time_steps * number of writers, if there is 1 group only
                */
                time_index[j][i] = (uint32_t *)
                    malloc (sizeof(uint32_t)*mh->pgs_count);
            } else {
                time_index[j][i] = (uint32_t *)
                    malloc (sizeof(uint32_t)*mh->time_steps);
            }
        }
    }

    root = &(fh->pgs_root);
    uint64_t grpid = grpidlist[0];
    uint32_t pg_time_count = 0, first_pg;
    uint32_t time_id = tidx_start;
    first_pg = 0; /* The first pg for a given timestep and group */
    for (i = 0; i < mh->pgs_count; i++) {
        pg_pids [i] = (*root)->process_id;
        pg_offsets [i] = (*root)->offset_in_file;
        //printf ("### root->time_index=%d,  time_id=%d\n", (*root)->time_index, time_id);
        if ((*root)->time_index == time_id) {
            /* processing still the same timestep */
            if (grpid == grpidlist[i]) {
                /* processing still the same group */
                /* FIXME: is this the order in the file? time..groups or group..times? */
                pg_time_count += 1;
            } else {
                /* changing group: pg_time_count is for the current group the number of pgs of the same time */
                time_index [0][grpid][time_id-tidx_start] = first_pg;
                time_index [1][grpid][time_id-tidx_start] = pg_time_count;
                //printf ("#-- time_index[0][%d][%d]=%d\n", grpid, time_id-tidx_start, first_pg);
                //printf ("#   time_index[1][%d][%d]=%d\n", grpid, time_id-tidx_start, pg_time_count);
                grpid = grpidlist [i];
                pg_time_count = 1;
                first_pg = i; // new group starts from this pg
            }
        }
        else {
            /* change in timestep */
            if (group_count == 1) {
                /* single group in file (the most frequent case) */
                time_index [0][grpid][time_id-tidx_start] = first_pg;
                time_index [1][grpid][time_id-tidx_start] = pg_time_count;
                //printf ("### time_index[0][%d][%d]=%d\n", grpid, time_id-tidx_start, first_pg);
                //printf ("    time_index[1][%d][%d]=%d\n", grpid, time_id-tidx_start, pg_time_count);
                first_pg = i;
            }
            else {
                if (grpid == grpidlist[i]) {
                    pg_time_count += 1;
                } else {
                    time_index [0][grpid][time_id-tidx_start] = first_pg;
                    time_index [1][grpid][time_id-tidx_start] = pg_time_count;
                    //printf ("#.. time_index[0][%d][%d]=%d\n", grpid, time_id-tidx_start, first_pg);
                    //printf ("    time_index[1][%d][%d]=%d\n", grpid, time_id-tidx_start, pg_time_count);
                    grpid = grpidlist [i];
                    first_pg = i;
                }
            }
            time_id = (*root)->time_index;
            pg_time_count = 1;
        }
        root = &(*root)->next;
    }
    /* set last grp/time count to complete the above procedure */
    time_index [0][grpid][time_id-tidx_start] = first_pg;
    time_index [1][grpid][time_id-tidx_start] = pg_time_count;
    //printf ("#   time_index[0][%d][%d]=%d\n", grpid, time_id-tidx_start, first_pg);
    //printf ("    time_index[1][%d][%d]=%d\n", grpid, time_id-tidx_start, pg_time_count);


    /* Copy group_count strings from namelist and then free up namelist */
    char ** grp_namelist;

    grp_namelist = (char **) malloc (sizeof(char*) * group_count);
    for (i=0;i<group_count;i++) {
        //grp_namelist[i] = (char *) malloc (strlen(namelist[i])+1);
        //strcpy(grp_namelist[i],namelist[i]);
        grp_namelist[i] = namelist[i];
    }
    free(namelist);

    // here we need:
    //        grp_namelist [ngroup]
    //    time_index   [2][ngroup][nprocess]
    //    pg_offsets   [npgs]

    free (pg_pids);

    fh->gvar_h = (struct BP_GROUP_VAR *) malloc (sizeof(struct BP_GROUP_VAR));
    fh->gvar_h->group_count = group_count;
    fh->gvar_h->pg_offsets = pg_offsets;
    fh->gvar_h->namelist = grp_namelist;
    fh->gvar_h->time_index = time_index;
    fh->gvar_h->group_id = 0;
    fh->gvar_h->var_offsets = 0;
    fh->gvar_h->var_namelist = 0;
    fh->gvar_h->var_counts_per_group = 0;

    fh->gattr_h = (struct BP_GROUP_ATTR *) malloc (sizeof(struct BP_GROUP_ATTR));
    fh->gattr_h->group_count = group_count;
    fh->gattr_h->namelist = grp_namelist;
    fh->gattr_h->group_id = 0;
    fh->gattr_h->attr_offsets = 0;
    fh->gattr_h->attr_namelist = 0;
    fh->gattr_h->attr_counts_per_group = 0;

    fh->tidx_start = tidx_start;
    fh->tidx_stop= tidx_stop;

    free(grpidlist);
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

