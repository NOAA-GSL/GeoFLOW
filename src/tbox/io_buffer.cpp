

#include "tbox/io_buffer.hpp"

namespace geoflow {
namespace tbox {


//
// Construct a I/O buffer object.  The object will require further
// initialization to set up I/O streams and the prefix string.
//
IOBuffer::IOBuffer(){
  setActive(true);
  setOutputStream1(NULL);
  setOutputStream2(NULL);
  setBufferSize(4096);
}

///
// The destructor deallocates internal data buffer.  It does not modify/
// the output streams.
//
IOBuffer::~IOBuffer(){
  if (d_buffer.size() > 0) {
    this->outputBuffer();
  }
  d_ostream1 = NULL;
  d_ostream2 = NULL;
}

//
// Activate or deactivate the output stream.
//
void IOBuffer::setActive(bool active){
  d_active = active;
}

//
// Set the prefix that begins every new line
//
void IOBuffer::setPrefixString(const std::string& text){
	d_prefix = text;
}

//
// Set the size of the buffer in Bytes
//
void IOBuffer::setBufferSize(const std::size_t& bytes){
	d_size = bytes;
	d_buffer.reserve(d_size);
}

//
// Set the primary output stream.
//
void IOBuffer::setOutputStream1(std::ostream* stream){
  d_ostream1 = stream;
}

//
// Set the secondary output stream.
//
void IOBuffer::setOutputStream2(std::ostream* stream){
  d_ostream2 = stream;
}

//
// Output a string to the output stream by invoking the
// outputString(string,length) method.
//
void IOBuffer::outputString(const std::string& text){
  outputString(text, text.size());
}

//
// Write a text string of the specified length to the output stream.
// Note that the string data is accumulated into the internal output
// buffer until an end-of-line is detected.
//
void IOBuffer::outputString(const std::string& text, const int length){
  if( (length > 0) && d_active ){

	  // Add prefix if needed
	  if( (d_buffer.length() == 0) && (d_prefix.length() > 0) ) {
		  copyToBuffer(d_prefix, d_prefix.size());
	  }

	  // If no new line character then copy entire text
	  // Else then recursively output lines
	  if( text.find('\n') == std::string::npos ){
		  copyToBuffer(text, length);

	  }
	  else {
		  const int ncopy = text.find('\n') + 1;
		  copyToBuffer(text, ncopy);
		  outputBuffer();
		  outputString(text.substr(ncopy), length - ncopy);
	  }

  }
}

//
// Copy data from the text string into the internal output buffer.
// If the internal buffer is not large enough to hold all of the string
// data, then fill and flush until complete.
//
void IOBuffer::copyToBuffer(const std::string& text, const int length){

  int num_copy;
  int start_ptr = 0;
  do {
    num_copy = d_size - d_buffer.length();
    
    d_buffer.append(text, start_ptr, num_copy);
    start_ptr = num_copy;

    if( start_ptr < length ){
      outputBuffer();
    }
  } while(start_ptr < length);

}

//
// Output buffered stream data to the active output streams and reset
// The buffer pointer to its empty state.
//
void IOBuffer::outputBuffer(){
  if( d_buffer.size() > 0 ){
    if (d_ostream1) {
      d_ostream1->write(d_buffer.c_str(),d_buffer.length());
    }
    if (d_ostream2) {
      d_ostream2->write(d_buffer.c_str(),d_buffer.length());
    }
    d_buffer.clear();
    d_buffer.reserve(d_size);
  }
}

//
// Synchronize the I/O buffer and write string data.  This routine
// is called from streambuf.
//
int IOBuffer::sync(){
  const int n = static_cast<int>(pptr() - pbase());
  if (n > 0) outputString(pbase(), n);
  return 0;
}

//
// Write the specified number of characters into the output stream.
// This routine is called from streambuf.  If this routine is not
// provided, then overflow() is called instead for each character.
//
// Note that this routine is not required; it only
// offers some efficiency over overflow().
//
#if !defined(__PGI) && !defined(__INTEL_COMPILER) && (defined(__GNUG__))
std::streamsize IOBuffer::xsputn(const std::string& text,
                                 std::streamsize n){
  sync();
  if (n > 0) outputString(text, static_cast<int>(n));
  return n;
}
#endif

//
// Write a single character into the I/O buffer.  This routine is
// called from streambuf.
//
int IOBuffer::overflow(int ch){
  const int n = static_cast<int>(pptr() - pbase());
  if (n && sync()) {
    return EOF;
  }
  if (ch != EOF) {
    char character[2];
    character[0] = (char)ch;
    character[1] = 0;
    outputString(character, 1);
  }
  pbump(-n);
  return 0;
}

} // namespace tbox
} // namespace geoflow
