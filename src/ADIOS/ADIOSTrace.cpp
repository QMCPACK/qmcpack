#ifdef HAVE_ADIOS
#include <cstdlib>
#include <boost/lexical_cast.hpp>
#include <Configuration.h>
#include "ADIOSTrace.h"


namespace ADIOS
{

    Trace::Trace (string group_name, MPI_Comm comm, xmlNodePtr adios_options)
    {
        group = 0;
        scalars_size = 0;
        rowsize = 0;
        transforms.clear();
        groupname = group_name;
        method_name = "MPI";
        method_args = "";
        if (adios_options)
            process_options (adios_options);

        adios_declare_group (&group, groupname.c_str(), "", adios_flag_yes);
        adios_select_method (group, method_name.c_str(), method_args.c_str(), "");
        // define some scalar vars 
        // first (global) dimension of all output arrays from all processes together
        adios_define_var (group, "nrows","", adios_integer, 0, 0, 0);
        scalars_size += sizeof(int);
        // first (local) dimension of all output arrays from one process
        adios_define_var (group, "/aux/max_nrows","", adios_integer, 0, 0, 0);
        scalars_size += sizeof(int);
        // first offset of a process in the global array
        adios_define_var (group, "/aux/offset","", adios_integer, 0, 0, 0);
        scalars_size += sizeof(int);
        // Information: first actual local dimension of all output arrays from one process
        // we fill all rows with 0 beyond this to read max_nrows_local rows on
        // all local processes
        adios_define_var (group, "/aux/actual_nrows_local","", adios_integer, 0, 0, 0);
        scalars_size += sizeof(int);

        adios_define_attribute (group, "nrows/description","", adios_string, 
            "Global dimension of the data arrays. Note that some rows are all-zero!", 0);
        adios_define_attribute (group, "/aux/max_nrows/description","", adios_string, 
            "Local (per MPI task) dimension of the data arrays. "
            "Equal on all tasks but some rows are all-zero!", 0);
        adios_define_attribute (group, "/aux/offset/description","", adios_string, 
            "Offset of each MPI task in the global array. It should be rank*max_nrows", 0);
        adios_define_attribute (group, "/aux/actual_nrows_local/description","", adios_string, 
            "The number of non-zero rows on each MPI task", 0);

        t_comm = comm;
        f = 0;
    }

    Trace::~Trace()
    {
        adios_free_group (group);
        transforms.clear();
    }

    int Trace::define_var(string path, int ndim, int *dims, string type)
    {
        int asize;
        enum ADIOS_DATATYPES atype;
        if (!type.compare("int"))
        {
            atype = adios_integer;    
            asize = sizeof (int);
        }
        else if (!type.compare("real"))
        {
            atype = adios_real;    
            asize = sizeof (float);
        }
        else if (!type.compare("double"))
        {
            atype = adios_double;    
            asize = sizeof (double);
            type = "double";
        }
        else if (!type.compare("complex"))
        {
            atype = adios_complex; 
            asize = sizeof (complex<float>);
            type = "complex";
        }
        else if (!type.compare("complex double"))
        {
            atype = adios_double_complex; 
            asize = sizeof (complex<double>);
            type = "complex double";
        }
        else 
        {
            cerr << "ADIOS:Trace:define_var: Wrong type is given:"<<type<< endl;
        }

        /* Scalars are not supported here */
        if (ndim < 1) 
        {
            cerr << "ADIOS:Trace:define_var: array is expected but ndim=1 for var "<<path<< endl;
            return 1;
        }

        string ldims = "/aux/max_nrows"; 
        string gdims = "nrows"; 
        string offs  = "/aux/offset"; 

        /* Single column arrays stay single column (1D array)
         * The rest will get the first dimension as extra dimension
         */
        if (ndim > 1 || dims[0] > 1) 
        {
            for (int i=0; i<ndim; i++) 
            {
                ldims = ldims + "," + boost::lexical_cast<std::string>(dims[i]); 
                gdims = gdims + "," + boost::lexical_cast<std::string>(dims[i]); 
                offs  = offs  + ",0";
                asize *= dims[i];
            }
        }

        qmcplusplus::app_log()<<"  ADIOS define: "<<type<<"\t"<<path<<"  ["<<ldims<<"]"<<endl;

        int64_t varid = adios_define_var (group, path.c_str(), "", atype, 
                                    ldims.c_str(), gdims.c_str(), offs.c_str());

        std::map<std::string,std::string>::iterator iter;
        iter = transforms.find(path);
        if (iter != transforms.end()) {
            qmcplusplus::app_log()<<"        set transform to: "<<iter->second.c_str()<<endl;
            adios_set_transform (varid, iter->second.c_str());
        }

        rowsize += asize;
    }

    void Trace::open (string filename)
    {
        //qmcplusplus::app_log()<<"  ADIOS open: "<<filename<<endl;
        int err = adios_open(&f, groupname.c_str(), filename.c_str(), "w", t_comm);
        if (err) {
            cerr << "ADIOS:Trace:open failed: "<<adios_get_last_errmsg()<< endl;
        }
    }

    void Trace::set_group_size (int nrows)
    {
        uint64_t group_size = scalars_size + (uint64_t)nrows * rowsize;
        uint64_t total_size;
        if (f) {
            //qmcplusplus::app_log()<<"  ADIOS group size: "<<group_size<<" bytes"<<endl;
            adios_group_size (f, group_size, &total_size);
            //qmcplusplus::app_log()<<"  ADIOS group total size: "<<total_size<<" bytes"<<endl;
        }
    }

    void Trace::write (string varname, void *data)
    {
        if (f) {
            //qmcplusplus::app_log()<<"  ADIOS write: "<<varname<<" data="<<data<<endl;
            adios_write (f, varname.c_str(), data);
        }
    }

    void Trace::close()
    {
        if (f) {
            //qmcplusplus::app_log()<<"  ADIOS close: "<<endl;
            adios_close (f);
            f = 0;
            MPI_Barrier (t_comm);
            /*
            int rank;
            MPI_Comm_rank (t_comm, &rank);
            cout<<"  ADIOS close finished on rank "<<rank<<endl;
            */
        }
    }

    void Trace::process_options (xmlNodePtr adios_options)
    {
        qmcplusplus::app_log()<<"  Process ADIOS options for traces"<<endl;
        xmlNodePtr element = adios_options->children;
        while(element!=NULL)
        {
            string name((const char*)element->name);
            if(name=="method")
            {
                string name = "MPI";
                OhmmsAttributeSet eattrib;
                eattrib.add(name,"name");
                eattrib.put(element);
                method_name = name;
                if (element->doc) {
                    const char *args = (const char*)
                        (xmlNodeListGetString(element->doc, element->xmlChildrenNode, 1));
                    if (args != NULL)  {
                        method_args = args;
                    }
                }
                qmcplusplus::app_log()<<"    Output method: "<<method_name<<endl;
                qmcplusplus::app_log()<<"        Arguments: "<<method_args<<endl;
            }
            else if(name=="var")
            {
                string name = "";
                string transform = "none";
                OhmmsAttributeSet eattrib;
                eattrib.add(name,"name");
                eattrib.add(transform,"transform");
                eattrib.put(element);
                if (name != "") {
                    transforms[name] = transform;
                    qmcplusplus::app_log()<<"    Transform "<<name<<"with "<<transform<<endl;
                } else {
                    qmcplusplus::app_log()<<"      No name given. Skip."<<endl;
                } 
            }
            else if (name!="text" && name!="comment")
            {
                APP_ABORT("ADIOS::Trace::process_options: "+name+" is not a valid sub-element of <adios_options> under <traces>\n  valid options are: method, var");
            }
            element=element->next;
        }

    }

    } //namespace

#endif // HAVE_ADIOS
