#ifndef OHMMS_QMC_WALKER_tf_H
#define OHMMS_QMC_WALKER_tf_H

herr_t file_info(hid_t loc_id,const char* name, void* opdata){
  H5G_stat_t statbuf;
  int* p = (int*)opdata;
  H5Gget_objinfo(loc_id,name,0,&statbuf);
  if(statbuf.type == H5G_GROUP) (*p)++;
}

#endif
