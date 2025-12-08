!******************************************************************************!
!                               zlib.f90
!------------------------------------------------------------------------------!
!> Fortran binding for the `zlib.h` header from the zlib file compression
!> library. No interface to gzprintf is included because there is no way to
!> call variadic functions from Fortran reliably. These interface are based on
!> the version 1.3.1 of the `zlib` library. The comments are from the `zlib.h`
!> header and modified to work with `FORD`. They are subject to the original
!> zlib licence:
!>
!>   Copyright (C) 1995-2024 Jean-loup Gailly and Mark Adler
!>
!>  This software is provided 'as-is', without any express or implied
!>  warranty.  In no event will the authors be held liable for any damages
!>  arising from the use of this software.
!>
!>  Permission is granted to anyone to use this software for any purpose,
!>  including commercial applications, and to alter it and redistribute it
!>  freely, subject to the following restrictions:
!>
!>  1. The origin of this software must not be misrepresented; you must not
!>     claim that you wrote the original software. If you use this software
!>     in a product, an acknowledgment in the product documentation would be
!>     appreciated but is not required.
!>  2. Altered source versions must be plainly marked as such, and must not be
!>     misrepresented as being the original software.
!>  3. This notice may not be removed or altered from any source distribution.
!>
!>  Jean-loup Gailly        Mark Adler
!>  jloup@gzip.org          madler@alumni.caltech.edu
!
!
!  Interfaces written by Paul Fossati, <paul.fossati@cea.fr>
!  Copyright (c) 2024-2025 CEA
!
!------------------------------------------------------------------------------!
! Redistribution and use in source and binary forms, with or without           !
! modification, are permitted provided that the following conditions are met:  !
!                                                                              !
!     * Redistributions of source code must retain the above copyright notice, !
!       this list of conditions and the following disclaimer.                  !
!                                                                              !
!     * Redistributions in binary form must reproduce the above copyright      !
!       notice, this list of conditions and the following disclaimer in the    !
!       documentation and/or other materials provided with the distribution.   !
!                                                                              !
!     * The name of the author may not be used to endorse or promote products  !
!      derived from this software without specific prior written permission    !
!      from the author.                                                        !
!                                                                              !
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"  !
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE    !
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE   !
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE     !
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR          !
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF         !
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS     !
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN      !
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)      !
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE   !
! POSSIBILITY OF SUCH DAMAGE.                                                  !
!******************************************************************************!
module zlib
    use iso_c_binding, only: gzFile => c_ptr, c_ptr, c_int, c_char, c_size_t
    implicit none(external, type)
    public

    ! Allowed flush values; see [[deflate(function)]] and [[inflate(function)]] for details.
    integer, parameter :: Z_NO_FLUSH      = 0
    integer, parameter :: Z_PARTIAL_FLUSH = 1
    integer, parameter :: Z_SYNC_FLUSH    = 2
    integer, parameter :: Z_FULL_FLUSH    = 3
    integer, parameter :: Z_FINISH        = 4
    integer, parameter :: Z_BLOCK         = 5
    integer, parameter :: Z_TREES         = 6

   ! Return codes for the compression/decompression functions. Negative values
   ! are errors, positive values are used for special but normal events.
   integer, parameter :: Z_OK            =  0
    integer, parameter :: Z_STREAM_END    =  1
    integer, parameter :: Z_NEED_DICT     =  2
    integer, parameter :: Z_ERRNO         = -1
    integer, parameter :: Z_STREAM_ERROR  = -2
    integer, parameter :: Z_DATA_ERROR    = -3
    integer, parameter :: Z_MEM_ERROR     = -4
    integer, parameter :: Z_BUF_ERROR     = -5
    integer, parameter :: Z_VERSION_ERROR = -6

    ! Compression levels
    integer, parameter :: Z_NO_COMPRESSION      =  0
    integer, parameter :: Z_BEST_SPEED          =  1
    integer, parameter :: Z_BEST_COMPRESSION    =  9
    integer, parameter :: Z_DEFAULT_COMPRESSION = -1

    ! compression strategy; see [[deflateInit2(function)]] for details
    integer, parameter :: Z_FILTERED         = 1
    integer, parameter :: Z_HUFFMAN_ONLY     = 2
    integer, parameter :: Z_RLE              = 3
    integer, parameter :: Z_FIXED            = 4
    integer, parameter :: Z_DEFAULT_STRATEGY = 0

    ! Possible values of the `data_type` field for [[deflate(function)]]
    integer, parameter :: Z_BINARY   = 0
    integer, parameter :: Z_TEXT     = 1
    integer, parameter :: Z_ASCII    = Z_TEXT
    integer, parameter :: Z_UNKNOWN  = 2

    ! The deflate compression method (the only one supported in this version)
    integer, parameter :: Z_DEFLATED = 8

    ! For initializing `zalloc`, `zfree`, `opaque`
    integer, parameter :: Z_NUL = 0

! Basic functions
    interface
!******************************************************************************!
!> The application can compare [[function(zlibVersion)]] and `ZLIB_VERSION`
!> for consistency.
!>
!> If the first character differs, the library code actually used is not
!> compatible with the `zlib.h` header file used by the application. This check
!> is automatically made by [[deflateInit(function)]] and
!> [[inflateInit(function)]].
!>
!>@note
!> The above comment comes from the C implementation. For now, `ZLIB_VERSION`
!> is not available from Fortran code.
!>@endnote
!******************************************************************************!
        type(c_ptr) function zlibVersion() bind(C, name="zlibVersion")
            import
            implicit none(external, type)
        end function
!******************************************************************************!
!> Initializes the internal stream state for compression.
!>
!> The fields `zalloc`, `zfree` and opaque must be initialized before by the
!> caller. If `zalloc` and `zfree` are set to `Z_NULL`,
!> [[deflateInit(function)]] updates them to use default allocation functions.
!> `total_in`, `total_out`, `adler`, and `msg` are initialized.
!>
!> The compression level must be `Z_DEFAULT_COMPRESSION`, or between `0` and
!> `9`: `1` gives best speed, `9` gives best compression, `0` gives no
!> compression at all (the input data is simply copied a block at a time).
!> `Z_DEFAULT_COMPRESSION` requests a default compromise between speed and
!> compression (currently equivalent to level `6`).
!>
!> [[deflateInit(function)]] returns `Z_OK` if success, `Z_MEM_ERROR` if
!> there was not enough memory, `Z_STREAM_ERROR` if `level` is not a valid
!> compression level, or `Z_VERSION_ERROR` if the `zlib` library version
!> ([[zlib_version(function)]]) is incompatible with the version assumed by the
!> caller (`ZLIB_VERSION`). `msg` is set to `null` if there is no error message.
!> [[deflateInit(function)]] does not perform any compression: this will be
!> done by [[deflate(function)]].
!>
!>@note
!> The `z_stream` derived type has no Fortran interface. At the moment, `strm`
!> is an opaque pointer.
!>@endnote
!******************************************************************************!
        integer(c_int) function deflateInit(strm, level) bind(C, name="deflateInit")
            import
            implicit none(external, type)
            type(c_ptr), value :: strm
            integer(c_int), value :: level
        end function

        integer(c_int) function deflate(strm, flush) bind(C, name="deflate")
            import
            implicit none(external, type)
            type(c_ptr), value :: strm
            integer(c_int), value :: flush
        end function
!******************************************************************************!
!> All dynamically allocated data structures for this stream are freed.
!>
!> This function discards any unprocessed input and does not flush any pending
!> output.
!>
!> [[deflateEnd(function)]] returns `Z_OK` if success, `Z_STREAM_ERROR` if the
!> stream state was inconsistent, `Z_DATA_ERROR` if the stream was freed
!> prematurely (some input or output was discarded). In the error case, `msg`
!> may be set but then points to a static string (which must not be
!> deallocated).
!>
!>@note
!> The `z_stream` derived type has no Fortran interface. At the moment, `strm`
!> is an opaque pointer.
!>@endnote
!******************************************************************************!
        integer(c_int) function deflateEnd(strm) bind(C, name="deflateEnd")
            import
            implicit none(external, type)
            type(c_ptr), value :: strm
        end function
!******************************************************************************!
!> Initializes the internal stream state for decompression.
!>
!> The fields `next_in`, `avail_in`, `zalloc`, `zfree` and `opaque` must be
!> initialized before by the caller. In the current version of
!> [[inflate(function)]], the provided input is not read or consumed. The
!> allocation of a sliding window will be deferred to the first call of
!> [[inflate(function)]] (if the decompression does not complete on the
!> first call). If `zalloc` and `zfree` are set to `Z_NULL`,
!> [[inflateInit(function)]] updates them to use default allocation functions.
!> `total_in`, `total_out`, `adler`, and `msg` are initialized.
!>
!> [[inflateInit(function)]] returns `Z_OK` if success, `Z_MEM_ERROR` if there
!> was not enough memory, `Z_VERSION_ERROR` if the zlib library version is
!> incompatible with the version assumed by the caller, or `Z_STREAM_ERROR` if
!> the parameters are invalid, such as a null pointer to the structure. `msg`
!> is set to `null` if there is no error message. [[inflateInit(function)]]
!> does not perform any decompression. Actual decompression will be done by
!> [[inflate(function)]]. So `next_in`, and `avail_in`, `next_out`, and
!> `avail_out` are unused and unchanged. The current implementation of
!> [[inflateInit(function)]] does not process any header information --
!> that is deferred until [[inflate(function)]] is called.
!>
!>@note
!> The `z_stream` derived type has no Fortran interface. At the moment, `strm`
!> is an opaque pointer.
!>@endnote
!******************************************************************************!
        integer(c_int) function inflateInit(strm) bind(C, name="inflateInit")
            import
            implicit none(external, type)
            type(c_ptr), value :: strm
        end function

        integer(c_int) function inflate(strm, flush) bind(C, name="inflate")
            import
            implicit none(external, type)
            type(c_ptr), value :: strm
            integer(c_int), value :: flush
        end function
!******************************************************************************!
!> All dynamically allocated data structures for this stream are freed.
!>
!> This function discards any unprocessed input and does not flush any pending
!> output.
!>
!> [[inflateEnd(function)]] returns `Z_OK` if success, or `Z_STREAM_ERROR` if
!> the stream state was inconsistent.
!>
!>@note
!> The `z_stream` derived type has no Fortran interface. At the moment, `strm`
!> is an opaque pointer.
!>@endnote
!******************************************************************************!
        integer(c_int) function inflateEnd(strm) bind(C, name="inflateEnd")
            import
            implicit none(external, type)
            type(c_ptr), value :: strm
        end function
    end interface

    ! File access functions
    interface
!******************************************************************************!
!> Open the gzip (.gz) file at path for reading and decompressing, or
!> compressing and writing.
!>
!> The mode parameter is as in `fopen` (`"rb"` or `"wb"`) but can also include
!> a compression level (`"wb9"`) or a strategy: `"f"` for filtered data as in
!> `"wb6f"`, `"h"` for Huffman-only compression as in `"wb1h"`, `"R"` for
!> run-length encoding as in `"wb1R"`, or `"F"` for fixed code compression
!> as in `"wb9F"` (See the description of [[deflateInit2(function)]] for more
!> information about the strategy parameter.). `"T"` will request transparent
!> writing or appending with no compression and not using the gzip format.
!>
!> `"a"` can be used instead of `"w"` to request that the gzip stream that will
!> be written be appended to the file. `"+"` will result in an error, since
!> reading and writing to the same gzip file is not supported. The addition of
!> `"x"` when writing will create the file exclusively, which fails if the file
!> already exists.  On systems that support it, the addition of "e" when
!> reading or writing will set the flag to close the file on an `execve()` call.
!>
!> These functions, as well as gzip, will read and decode a sequence of gzip
!> streams in a file. The append function of [[gzopen(function)]] can be used
!> to create such a file. (Also see [[gzflush(function)]] for another way to do
!> this). When appending, [[gzopen(function)]] does not test whether the file
!> begins with a gzip stream, nor does it look for the end of the gzip streams
!> to begin appending. [[gzopen(function)]] will simply append a gzip stream
!> to the existing file.
!>
!> [[gzopen(function)]] can be used to read a file which is not in gzip format;
!> in this case [[gzread(function)]] will directly read from the file without
!> decompression. When reading, this will be detected automatically by looking
!> for the magic two-byte gzip header.
!>
!> [[gzopen(function)]] returns `NULL` if the file could not be opened, if there
!> was insufficient memory to allocate the `gzFile` state, or if an invalid mode
!> was specified (an `"r"`, `"w"`, or `"a"` was not provided, or `"+"` was
!> provided). [[errno(function)]] can be checked to determine if the reason
!> [[gzopen(function)]] failed was that the file could not be opened.
!******************************************************************************!
        type(gzFile) function gzopen(path, mode) bind(C, name="gzopen")
            import
            implicit none(external, type)
            character(1,c_char), dimension(*), intent(in) :: path
            character(1,c_char), dimension(*), intent(in) :: mode
        end function
!******************************************************************************!
!> Associate a gzFile with the file descriptor `fd`.
!>
!> File descriptors are obtained from calls like `open`, `dup`, `creat`, `pipe`
!> or `fileno` (if the file has been previously opened with `fopen`). The mode
!> parameter is as in [[gzopen(function)]].
!>
!> The next call of [[gzclose(function)]] on the returned `gzFile` will also
!> close the file descriptor `fd`, just like `fclose(fdopen(fd, mode))` closes
!> the file descriptor `fd`.  If you want to keep fd open, use fd = dup(fd_keep); gz = gzdopen(fd,
!> mode);.  The duplicated descriptor should be saved to avoid a leak, since
!> gzdopen does not close fd if it fails.  If you are using fileno() to get the
!> file descriptor from a FILE *, then you will have to use dup() to avoid
!> double-close()ing the file descriptor.  Both gzclose() and fclose() will
!> close the associated file descriptor, so they need to have different file
!> descriptors.
!>
!> gzdopen returns NULL if there was insufficient memory to allocate the
!> gzFile state, if an invalid mode was specified (an 'r', 'w', or 'a' was not
!> provided, or '+' was provided), or if fd is -1.  The file descriptor is not
!> used until the next gz* read, write, seek, or close operation, so gzdopen
!> will not detect if fd is invalid (unless fd is -1).
!******************************************************************************!
        type(gzFile) function gzdopen(fd, mode) bind(C, name="gzdopen")
            import
            implicit none(external, type)
            integer(c_int), value :: fd
            character(1,c_char), dimension(*), intent(in) :: mode
        end function
!******************************************************************************!
!> Flush all pending output for `file`, if necessary, close `file` and
!> deallocate the (de)compression state.
!>
!> Note that once `file` is closed, you cannot call [[gzerror(function)]] with
!> `file`, since its structures have been deallocated. [[gzclose(function)]]
!> must not be called more than once on the same file, just as `free` must not
!> be called more than once on the same allocation.
!>
!> [[gzclose(function)]] will return `Z_STREAM_ERROR` if `file` is not valid,
!> `Z_ERRNO` on a file operation error, `Z_MEM_ERROR` if out of memory,
!> `Z_BUF_ERROR` if the last read ended in the middle of a gzip stream, or
!> `Z_OK` on success.
!******************************************************************************!
        integer(c_int) function gzclose(file) bind(C, name='gzclose')
            import
            implicit none(external, type)
            type(gzFile), value :: file
        end function
!******************************************************************************!
!> Return the error message for the last error which occurred on file.
!>
!> `errnum` is set to zlib error number. If an error occurred in the file system
!> and not in the compression library, `errnum` is set to `Z_ERRNO` and the
!> application may consult [[errno(function]]) to get the exact error code.
!>
!> The application must not modify the returned string. Future calls to
!> this function may invalidate the previously returned string. If `file` is
!> closed, then the string previously returned by gzerror will no longer be
!> available.
!>
!> `gzerror()` should be used to distinguish errors from end-of-file for those
!> functions above that do not distinguish those cases in their return values.
!******************************************************************************!
        type(c_ptr) function gzerror(file, errnum) bind(C, name='gzerror')
            import
            implicit none(external, type)
            type(gzFile), value :: file
            integer(c_int) :: errnum
        end function
!******************************************************************************!
!> Clear the error and end-of-file flags for `file`.
!>
!> This is analogous to the `clearerr()` function in `stdio`. This is useful for
!> continuing to read a gzip file that is being written concurrently.
!******************************************************************************!
        subroutine gzclearerr(file) bind(C, name="gzclearerr")
            import
            implicit none(external, type)
            type(gzFile), value :: file
        end subroutine
!******************************************************************************!
!> Return `1` if the end-of-file indicator for file has been set while reading,
!> `0` otherwise.
!>
!> Note that the end-of-file indicator is set only if the read tried to go past
!> the end of the input, but came up short. Therefore, just like `feof()`,
!> `gzeof()` may return false even if there is no more data to read, in the
!> event that the last read request was for the exact number of bytes remaining
!> in the input file. This will happen if the input file size is an exact
!> multiple of the buffer size.
!>
!> If [[gzeof(function)]] returns `1`, then the read functions will return no
!> more data, unless the end-of-file indicator is reset by
!> [[gzclearerr(function)]] and the input file has grown since the previous end
!> of file was detected.
!******************************************************************************!
        integer(c_int) function gzeof(file) bind(C, name='gzeof')
            import
            implicit none(external, type)
            type(gzFile),   value :: file
        end function
!******************************************************************************!
!> Set the internal buffer size used by this library's functions for file to
!> size.
!>
!> The default buffer size is 8192 bytes. This function must be called after
!> [[gzopen(function)]] or [[gzdopen(function)]], and before any other calls
!> that read or write the file. The buffer memory allocation is always deferred
!> to the first read or write. Three times that size in buffer space is
!> allocated. A larger buffer size of, for example, 64K or 128K bytes will
!> noticeably increase the speed of decompression (reading).
!>
!> The new buffer size also affects the maximum length for
!> [[gzprintf(function)]].
!>
!> [[gzbuffer(function)]] returns `0` on success, or `-1` on failure, such as
!> being called too late.
!******************************************************************************!
        integer(c_int) function gzbuffer(file, size) bind(C, name="gzbuffer")
            import
            implicit none(external, type)
            type(c_ptr), value :: file
            integer(c_int), value :: size !< Unsigned in the C library
        end function
!******************************************************************************!
!> Dynamically update the compression level and strategy for file.
!>
!> See the description of [[deflateInit2(function)]] for the meaning of these
!> parameters. Previously provided data is flushed before applying the parameter
!> changes.
!>
!> [[gzsetparams(function)]] returns `Z_OK` if success, `Z_STREAM_ERROR` if the
!> file was not opened for writing, `Z_ERRNO` if there is an error writing the
!> flushed data, or `Z_MEM_ERROR` if there is a memory allocation error.
!******************************************************************************!
        integer(c_int) function gzsetparams(file, level, strategy) bind(C, name="gzsetparams")
            import
            implicit none(external, type)
            type(gzFile), value :: file
            integer(c_int), value :: level
            integer(c_int), value :: strategy
        end function

        integer(c_int) function gzread(file, buff, len) bind(C, name="gzread")
            import
            implicit none(external, type)
            type(gzFile), value :: file
            character(1,c_char), dimension(*), intent(out) :: buff
            integer(c_int), value :: len !< Unsigned in the C library
        end function

        integer(c_size_t) function gzfread(buff, size, nitems, file) bind(C, name="gzfread")
            import
            implicit none(external, type)
            character(1,c_char), dimension(*), intent(out) :: buff
            integer(c_size_t), value :: size
            integer(c_size_t), value :: nitems
            type(gzFile), value :: file
        end function
!******************************************************************************!
!> Compress and write the len uncompressed bytes at `buf` to `file`.
!>
!> [[gzwrite(function)]] returns the number of uncompressed bytes written or
!> `0` in case of error.
!******************************************************************************!
        integer(c_int) function gzwrite(file, buf, len) bind(C, name="gzwrite")
            import
            implicit none(external, type)
            type(gzFile), value :: file
            character(1,c_char), dimension(*), intent(in) :: buf
            integer(c_int), value :: len !< Unsigned in the C library
        end function
!******************************************************************************!
!> Compress and write `nitems` items of size `size` from `buf` to `file`,
!> duplicating the interface of stdio's fwrite(), with `size_t` request and
!> return types.
!>
!> If the library defines `size_t`, then `z_size_t` is identical to `size_t`.
!> If not, then `z_size_t` is an unsigned integer type that can contain a
!> pointer.
!>
!> [[gzfwrite(function)]] returns the number of full items written of size
!> `size`, or zero if there was an error. If the multiplication of `size` and
!> `nitems` overflows, i.e. the product does not fit in a `z_size_t`, then
!> nothing is written, zero is returned, and the error state is set to
!> `Z_STREAM_ERROR`.
!******************************************************************************!
        integer(c_size_t) function gzfwrite(buf, size, nitems, file) bind(C, name="gzfwrite")
            import
            implicit none(external, type)
            character(1,c_char), dimension(*), intent(in) :: buf
            integer(c_size_t), value :: size
            integer(c_size_t), value :: nitems
            type(gzFile), value :: file
        end function
!******************************************************************************!
!> Compress and write the given null-terminated string `s` to `file`, excluding
!> the terminating null character.
!>
!> [[gzputs(function)]] returns the number of characters written, or `-1` in
!> case of error.
!******************************************************************************!
        integer(c_int) function gzputs(file, s) bind(C, name="gzputs")
            import
            implicit none(external, type)
            type(gzFile), value :: file
            character(1,c_char), dimension(*), intent(in) :: s
        end function
!******************************************************************************!
!> Read and decompress bytes from `file` into `buf`, until `len-1` characters
!> are read, or until a newline character is read and transferred to `buf`, or
!> an end-of-file condition is encountered.
!>
!> If any characters are read or if `len` is one, the string is terminated with
!> a null character. If no characters are read due to an end-of-file or `len`
!> is less than one, then the buffer is left untouched.
!>
!> [[gzgets(function)]] returns `buf` which is a null-terminated string, or it
!> returns `NULL` for end-of-file or in case of error. If there was an error,
!> the contents at buf are indeterminate.
!******************************************************************************!
        integer(c_int) function gzgets(file, buf, len) bind(C, name="gzgets")
            import
            implicit none(external, type)
            type(gzFile), value :: file
            character(1,c_char), dimension(*), intent(out) :: buf
            integer(c_int), value :: len !< Unsigned in the C library
        end function
!******************************************************************************!
!> Compress and write `c`, converted to an unsigned char, into `file`.
!>
!> [[gzputc(function)]] returns the value that was written, or `-1` in case of
!> error.
!******************************************************************************!
        integer(c_int) function gzputc(file, c) bind(C, name="gzputc")
            import
            implicit none(external, type)
            type(gzFile), value :: file
            integer(c_int), value :: c
        end function
!******************************************************************************!
!> Read and decompress one byte from `file`.
!>
!> [[gzgetc(function)]] returns this byte or `-1` in case of end of file or
!> error. This is implemented as a macro for speed. As such, it does not do all
!> of the checking the other functions do. I.e. it does not check to see if
!> `file` is `NULL`, nor whether the structure file points to has been
!> clobbered or not.
!******************************************************************************!
        integer(c_int) function gzgetc(file) bind(C, name="gzgetc")
            import
            implicit none(external, type)
            type(gzFile), value :: file
        end function
!******************************************************************************!
!> Push `c` back onto the stream for `file` to be read as the first character
!> on the next read.
!>
!> At least one character of push-back is always allowed.
!> [[gzungetc(function)]] returns the character pushed, or `-1` on failure.
!> [[gzungetc(function)]] will fail if `c` is `-1`, and may fail if a character
!> has been pushed but not read yet. If [[gzungetc(function)]] is used
!> immediately after [[gzopen(function)]] or [[gzdopen(function)]], at least the
!> output buffer size of pushed characters is allowed. (See
!> [[gzbuffer(function)]].) The pushed character will be discarded if the
!> stream is repositioned with [[gzseek(function)]] or [[gzrewind(function)]].
!******************************************************************************!
        integer(c_int) function gzungetc(c, file) bind(C, name="gzungetc")
            import
            implicit none(external, type)
            integer(c_int), value :: c
            type(gzFile), value :: file
        end function
!******************************************************************************!
!> Flush all pending output to file.
!>
!> The parameter flush is as in the [[deflate(function)]] function. The return
!> value is the zlib error number (see function [[gzerror(function)]]).
!> [[gzflush(function)]] is only permitted when writing.
!>
!> If the flush parameter is `Z_FINISH`, the remaining data is written and the
!> gzip stream is completed in the output. If [[gzwrite(function)]] is called
!> again, a new gzip stream will be started in the output. [[gzread(function)]]
!> is able to read such concatenated gzip streams.
!>
!> [[gzflush(function)]] should be called only when strictly necessary because
!> it will degrade compression if called too often.
!******************************************************************************!
        integer(c_int) function gzflush(file, flush) bind(C, name="gzflush")
            import
            implicit none(external, type)
            type(gzFile), value :: file
            integer(c_int), value :: flush
        end function
!******************************************************************************!
!> Set the starting position to offset relative to whence for the next
!> [[gzread(function)]] or [[gzwrite(function)]] on `file`.
!>
!> The offset represents a number of bytes in the uncompressed data stream.
!> The `whence` argument is defined as in `lseek`; the value `SEEK_END`
!> is not supported.
!>
!> If the file is opened for reading, this function is emulated but can be
!> extremely slow. If the file is opened for writing, only forward seeks are
!> supported; [[gzseek(function)]] then compresses a sequence of zeroes up to
!> the new starting position.
!>
!> [[gzseek(function)]] returns the resulting offset location as measured in
!> bytes from the beginning of the uncompressed stream, or `-1` in case of error,
!> in particular if the file is opened for writing and the new starting position
!> would be before the current position.
!******************************************************************************!
        integer(c_size_t) function gzseek(file, offset, whence) bind(C, name='gzseek')
            import
            implicit none(external, type)
             type(gzFile), value :: file
            integer(c_size_t), value :: offset
            integer(c_int), value :: whence
        end function
!******************************************************************************!
!> Rewind file.
!>
!> This function is supported only for reading.
!>
!> `gzrewind(file)` is equivalent to `gzseek(file, 0_c_long, SEEK_SET)`.
!******************************************************************************!
        integer(c_int) function gzrewind(file) bind(C, name='gzrewind')
            import
            implicit none(external, type)
            type(gzFile), value :: file
        end function
!******************************************************************************!
!> Return the starting position for the next [[gzread(function)]] or
!> [[gzwrite(function)]] on `file`.
!>
!> This position represents a number of bytes in the uncompressed data stream,
!> and is zero when starting, even if appending or reading a gzip stream from
!> the middle of a file using [[gzdopen(function)]].
!>
!> `gztell(file)` is equivalent to `gzseek(file, 0_c_long, SEEK_CUR)`
!******************************************************************************!
        integer(c_size_t) function gztell(file) bind(C, name='gztell')
            import
            implicit none(external, type)
            type(gzFile), value :: file
        end function
    end interface
end module
