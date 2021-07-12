/**
 * @file CustomExceptions.h
 * @author Till Budde (tilljanis.budde@tu-dortmund.de)
 * @brief Defines a number of custom exceptions using BOOST
 * @version 0.1
 * @date 2021-06-28
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#ifndef CUSTOMEXCEPTIONS_H
#define CUSTOMEXCEPTIONS_H

#include<boost/exception/all.hpp>
#include<string>

namespace cexc {

  typedef boost::error_info<struct tag_errno_code,int> errno_code;
  /**
   * @brief Basic exceptoin
   * 
   */
  struct exception_base: virtual std::exception, virtual boost::exception { };
  /**
   * @brief Exception concerning in-/output
   * 
   */
  struct io_error: virtual exception_base { };
  /**
   * @brief Not implemented error
   * 
   */
  struct not_imp_error: virtual exception_base { };
  /**
   * @brief Material Model not implemented error
   * 
   */
  struct not_mat_error: virtual not_imp_error { };

  struct runtime_error: virtual exception_base{ };

  struct convergence_error: virtual runtime_error{};
  
  /**
   * @brief Error concerning files
   * 
   */
  struct file_error: virtual io_error { };
  /**
   * @brief File can not be read
   * 
   */
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
