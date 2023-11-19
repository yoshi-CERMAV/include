//
//  Read_D5WOS.h
//  D2AM2023
//
//  Created by yoshi on 13/09/2023.
//

#ifndef Read_D5WOS_h
#define Read_D5WOS_h
#include <hdf5.h>
#include <hdf5_hl.h>
#include <H5Cpp.h>
#include <cstring>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
using namespace std;
using namespace H5;
#define SIZE_WOS 600 * 1156
#define SIZE_D5 578 * 960
#define YN_MAX_LEN 256
class read_h5
{
public:
    read_h5(char filename[]){
        strncpy(filename_,filename, YN_MAX_LEN);
        strncpy(D5fmt_,"/%d.1/measurement/D5", YN_MAX_LEN);
        strncpy(WOSfmt_,"/%d.1/measurement/WOS", YN_MAX_LEN);
        file_ = new H5File(filename, H5F_ACC_RDONLY);
        file_set = true;
    }
    read_h5(){
        strncpy(path_,"/Users/yoshi/temp/putaux/",YN_MAX_LEN);
        strncpy(filename_, "sample3/sample3_0001/sample3_0001.h5",YN_MAX_LEN);
        strncpy(D5fmt_,"/%d.1/measurement/D5", YN_MAX_LEN);
        strncpy(WOSfmt_,"/%d.1/measurement/WOS", YN_MAX_LEN);
        cout << "D5fmt"<<D5fmt_<<endl;
//        char full_path[255];
//        strncpy(full_path, path_, 255);
//        strncat(full_path, filename_, 255);
//        file_ = new H5File(full_path, H5F_ACC_RDONLY);
    }

    int set_path(const char path[])
    {
        struct stat sb;
        if(stat(path, &sb)) {cerr<<"directory does not exist"; return 0;}
        else strncpy(path_, path, YN_MAX_LEN);
        return 0;
    }
    int set_filename(const char filename[])
    {
  //      char full_path[255];
        strncpy(filename_,filename, YN_MAX_LEN);
  //      strncpy(full_path, path_, 255);
  //      strncat(full_path, filename_, 255);
        struct stat buf;
        if(stat(filename, &buf)) {cerr<<"directory does not exist"; return 1;}
        if(!stat(filename, &buf)){
            if(file_) delete file_;
            file_ = new H5File(filename, H5F_ACC_RDONLY);
       }
        cout <<"fmt"<<D5fmt_<<endl;
        file_set = true;
        return 0;
    }
    
    int read_D5(int i, int *data){
        cout << "reading D5"<<endl;
        return read_data(i, data, D5fmt_, SIZE_D5*10);
    }
    int read_WOS(int i, int *data){
        return read_data(i, data, WOSfmt_, SIZE_WOS*10);
    }
    int read_data(int i, int *data, char fmt[], size_t max_size){
        cout << "reading"<<endl;
        if(file_set==false){
            cerr<< "file not yet set";
            cout << filename_;return 1;
        }
        char name [YN_MAX_LEN];
        cout <<"fmt"<< fmt<<" "<<D5fmt_<<endl;
        snprintf(name, YN_MAX_LEN, fmt, i);
        cout <<"name "<<name <<endl;
        DataSet dataset;
        cout << file_<<endl;
        struct stat buf;

        if (!file_){
            if(!stat(filename_, &buf)){
                file_ = new H5File(filename_, H5F_ACC_RDONLY);
                
            }else cerr << filename_ <<"does not exist"<<endl;
        }
        try {dataset = file_->openDataSet(name);}
        catch(FileIException &E){cerr <<"could not open dataspace";return 2;}
        if( !dataset.getId() )
            {
                cerr<<"ReportReaderHDF5: "
                     <<"Dataset " << name << " not found "
                     <<"in file: " << file_->getFileName();
                return -1;
            }
        cout <<"get Space"<<endl;
        DataSpace dataspace = dataset.getSpace();
        cout <<"space OK"<<endl;
        size_t required_size = dataset.getStorageSize();
        int rank = dataspace.getSimpleExtentNdims();
        cout <<"rank "<< rank;
        required_size =dataspace.getSimpleExtentNpoints();
        cout << "size "<<required_size<<endl;;
        cout << "required_size "<< required_size<<" "<< max_size<<endl;
        cout << required_size / max_size <<endl;
        if(required_size > max_size){cerr <<"not enough space "<<required_size  ; return 0;}
        if(!data) {cout<< "no memory "<<endl; return 1;}
        dataset.read(data, PredType::STD_I32LE, dataspace, dataspace);
        return 0;
    }
    ~read_h5(){file_->close();}
    char * filename(){return filename_;}
private:
    char path_[YN_MAX_LEN];
    char filename_[YN_MAX_LEN];
    char D5fmt_[YN_MAX_LEN];
    char WOSfmt_[YN_MAX_LEN];
    bool file_set = false;
    H5File *file_;
};

#endif /* Read_D5WOS_h */
