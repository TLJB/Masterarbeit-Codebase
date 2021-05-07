#ifndef CUSTOMEXCEPTIONS_H
#define CUSTOMEXCEPTIONS_H

#include<boost/exception/all.hpp>
#include<string>

namespace cexc {

  typedef boost::error_info<struct tag_errno_code,int> errno_code;
  struct exception_base: virtual std::exception, virtual boost::exception { };
  struct io_error: virtual exception_base { };
  struct file_error: virtual io_error { };
  struct read_error: virtual io_error { };
  struct file_read_error: virtual file_error, virtual read_error {
    public:
      file_read_error(){
        message = "Default message";
      }
      file_read_error(std::string msg){
        message = msg;
      }
      std::string message;
      const char* what() const throw(){
        return message.c_str();
      }

   };

}


#endif 
