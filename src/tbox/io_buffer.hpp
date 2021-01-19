
#ifndef IOBUFFER_H_
#define IOBUFFER_H_


#include <iostream>
#include <string>

namespace geoflow {
namespace tbox {


/**
 * Class IOBuffer is a simple I/O stream utility that
 * intercepts output from an ostream and redirects the output
 * for I/O.  This class defines a stream buffer class for an
 * ostream class.
 */
class IOBuffer : public std::streambuf{
public:

   /**
    * Create a I/O buffer class.  The object will require further
    * initialization to set up the I/O streams and prefix string.
    */
   IOBuffer();

   /**
    * The destructor simply deallocates any internal data
    * buffers.  It does not modify the output streams.
    */
   virtual ~IOBuffer();

   /**
    * Set whether the output stream will be active.  If the I/O buffer
    * stream is disabled, then no data is forwarded to the output streams.
    * The internal data buffer is deallocated and pointers are reset
    * whenever the I/O buffer is deactivated.
    */
   void setActive(bool active);

   /**
    * Set the maximum size of the internal buffer for storing characters.
    * When the size of the internal buffer is reached a flush occurs to the
    * streams clearing the way for more data to be stored.
    */
   void setBufferSize(const std::size_t& bytes);

   /**
    * Set the prefix that begins every new line to the output stream.
    * A sample prefix is "P=XXXXX: ", where XXXXX represents the MPI
    * rank.
    */
   void setPrefixString(const std::string& text);

   /**
    * Set the primary output stream.  If not NULL, then output data is
    * sent to this stream.  The primary output stream is typically stderr
    * or stdout or perhaps a log file.
    */
   void setOutputStream1(std::ostream* stream);

   /**
    * Set the secondary output stream.  If not NULL, then output data is sent
    * to this stream.  The secondary output stream is typically NULL or a log
    * file that mirrors the primary output stream.
    */
   void setOutputStream2(std::ostream* stream);

   /**
    * Write a text string to the output stream.  Note that the string is
    * not actually written until an end-of-line is detected.
    */
   void outputString(const std::string& text);

   /**
    * Write a text string of the specified length to the output file.  Note
    * that the string is not actually written until an end-of-line is detected.
    */
   void outputString(const std::string& text, const int length);

   /**
    * Synchronize the I/O buffer (called from streambuf).
    */
   virtual int sync();

#if !defined(__PGI) && !defined(__INTEL_COMPILER) && (defined(__GNUG__))
   /**
    * Write the specified number of characters into the output stream (called
    * from streambuf).
    */
   virtual std::streamsize xsputn(const std::string& text,std::streamsize n);
#endif

   /**
    * Write an overflow character into the I/O buffer (called from
    * streambuf).
    */
   virtual int overflow(int ch);


private:
   void copyToBuffer(const std::string& text,const int length);
   void outputBuffer();          // output internal buffer data to streams

   std::ostream* d_ostream1;     // primary output stream for buffer
   std::ostream* d_ostream2;     // secondary output stream (e.g., for log file)
   std::string   d_buffer;       // internal buffer to store accumulated string
   std::string   d_prefix;       // prefix string added to start of each line
   std::size_t   d_size;         // max buffer size in Bytes
   bool d_active;                // whether this output stream is active
};

} // namespace tbox
} // namespace geoflow

#endif
