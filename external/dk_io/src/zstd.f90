!******************************************************************************!
!                              zstd.f90
!------------------------------------------------------------------------------!
!> Fortran binding for functions in the `zstd.h` header from the zstandard file
!> compression library.
!
!  Written by Paul Fossati, <paul.fossati@cea.fr>
!  Copyright (c) 2024-2025 CEA
!
!------------------------------------------------------------------------------!
! Redistribution and use in source and binary forms, with or without           !
! modification, are permitted provided that the following conditions are met:  !
!                                                                              !
!    !> Redistributions of source code must retain the above copyright notice, !
!       this list of conditions and the following disclaimer.                  !
!                                                                              !
!    !> Redistributions in binary form must reproduce the above copyright      !
!       notice, this list of conditions and the following disclaimer in the    !
!       documentation and/or other materials provided with the distribution.   !
!                                                                              !
!    !> The name of the author may not be used to endorse or promote products  !
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
module zstd
    use iso_c_binding
    use iso_c_binding, only: ZSTD_CCtx => c_ptr, ZSTD_DCtx => c_ptr, &
        ZSTD_CStream => c_ptr, ZSTD_DStream => c_ptr, c_ptr

    implicit none(external, type)
    public

    ! Default constant
    integer, parameter :: ZSTD_CLEVEL_DEFAULT = 3

    ! Constants
    integer, parameter :: ZSTD_MAGICNUMBER = int(z'0FD2FB528')
    integer, parameter :: ZSTD_MAGIC_DICTIONARY = int(z'EC30A437')
    integer, parameter :: ZSTD_MAGIC_SKIPPABLE_START = int(z'184D2A50') !< all 16 values, from 0x184D2A50 to 0x184D2A5F, signal the beginning of a skippable frame
    integer, parameter :: ZSTD_MAGIC_SKIPPABLE_MASK = int(z'FFFFFFF0')

    integer, parameter :: ZSTD_BLOCKSIZELOG_MAX = 17
    integer, parameter :: ZSTD_BLOCKSIZE_MAX = ishft(1, ZSTD_BLOCKSIZELOG_MAX)

    enum, bind(C)
        enumerator :: ZSTD_fast     = 1
        enumerator :: ZSTD_dfast    = 2
        enumerator :: ZSTD_greedy   = 3
        enumerator :: ZSTD_lazy     = 4
        enumerator :: ZSTD_lazy2    = 5
        enumerator :: ZSTD_btlazy2  = 6
        enumerator :: ZSTD_btopt    = 7
        enumerator :: ZSTD_btultra  = 8
        enumerator :: ZSTD_btultra2 = 9
    end enum

    ! Version
    interface
!******************************************************************************!
!> Return runtime library version, the value is
!> `(MAJOR*100*100 + MINOR*100 + RELEASE)`.
!******************************************************************************!
        integer(c_int) function ZSTD_versionNumber() bind(C, name="ZSTD_versionNumber")
            import
            implicit none(external, type)
        end function
!******************************************************************************!
!> Return runtime library version, like "1.4.5". Requires v1.3.0+.
!******************************************************************************!
        type(c_ptr) function ZSTD_versionString() bind(C, name="ZSTD_versionString")
            import
            implicit none(external, type)
        end function
    end interface

    ! Simple Core API
    interface
!******************************************************************************!
!> Compresses `src` content as a single zstd compressed frame into already
!> allocated `dst`.
!>
!>@note
!> Providing `dstCapacity >= ZSTD_compressBound(srcSize)` guarantees that zstd
!> will have enough space to successfully compress the data.
!>@endnote
!>
!>@return : compressed size written into `dst` (<= `dstCapacity), or an error
!> code if it fails (which can be tested using [[ZSTD_isError(function)]]).
!******************************************************************************!
        integer(c_size_t) function ZSTD_compress(dst, dstCapacity, src, srcSize, compressionLevel) &
            bind(C, name="ZSTD_compress")

            use iso_c_binding
            character(1,c_char), dimension(*), intent(out) :: dst
            integer(c_size_t), value :: dstCapacity
            character(1,c_char), dimension(*), intent(in) :: src
            integer(c_size_t), value :: srcSize
            integer(c_int), value :: compressionLevel
        end function
!******************************************************************************!
!> `compressedSize` : must be the _exact_ size of some number of compressed and/or skippable frames.
!>  Multiple compressed frames can be decompressed at once with this method.
!>  The result will be the concatenation of all decompressed frames, back to back.
!> `dstCapacity` is an upper bound of originalSize to regenerate.
!>  First frame's decompressed size can be extracted using ZSTD_getFrameContentSize().
!>  If maximum upper bound isn't known, prefer using streaming mode to decompress data.
!> @return : the number of bytes decompressed into `dst` (<= `dstCapacity`),
!>           or an errorCode if it fails (which can be tested using ZSTD_isError()). */
!******************************************************************************!
        integer(c_size_t) function ZSTD_decompress(dst, dstCapacity, src, compressedSize) bind(C, name="ZSTD_decompress")
            use iso_c_binding
            type(c_ptr), value :: dst
            integer(c_size_t), value :: dstCapacity
            type(c_ptr), value :: src
            integer(c_size_t), value :: compressedSize
        end function
!******************************************************************************!
!> `src` should point to the start of a ZSTD encoded frame.
!> `srcSize` must be at least as large as the frame header.
!>           hint : any size >= `ZSTD_frameHeaderSize_max` is large enough.
!> @return : - decompressed size of `src` frame content, if known
!>           - ZSTD_CONTENTSIZE_UNKNOWN if the size cannot be determined
!>           - ZSTD_CONTENTSIZE_ERROR if an error occurred (e.g. invalid magic number, srcSize too small)
!>  note 1 : a 0 return value means the frame is valid but "empty".
!>           When invoking this method on a skippable frame, it will return 0.
!>  note 2 : decompressed size is an optional field, it may not be present (typically in streaming mode).
!>           When `return==ZSTD_CONTENTSIZE_UNKNOWN`, data to decompress could be any size.
!>           In which case, it's necessary to use streaming mode to decompress data.
!>           Optionally, application can rely on some implicit limit,
!>           as ZSTD_decompress() only needs an upper bound of decompressed size.
!>           (For example, data could be necessarily cut into blocks <= 16 KB).
!>  note 3 : decompressed size is always present when compression is completed using single-pass functions,
!>           such as ZSTD_compress(), ZSTD_compressCCtx() ZSTD_compress_usingDict() or ZSTD_compress_usingCDict().
!>  note 4 : decompressed size can be very large (64-bits value),
!>           potentially larger than what local system can handle as a single memory segment.
!>           In which case, it's necessary to use streaming mode to decompress data.
!>  note 5 : If source is untrusted, decompressed size could be wrong or intentionally modified.
!>           Always ensure return value fits within application's authorized limits.
!>           Each application can set its own limits.
!>  note 6 : This function replaces ZSTD_getDecompressedSize() */
!******************************************************************************!
        integer(c_long_long) function ZSTD_getFrameContentSize(src, srcSize) bind(C, name="ZSTD_getFrameContentSize")
            use iso_c_binding
            type(c_ptr), value :: src
            integer(c_size_t), value :: srcSize
        end function
!******************************************************************************!
!>  This function is now obsolete, in favor of ZSTD_getFrameContentSize().
!>  Both functions work the same way, but ZSTD_getDecompressedSize() blends
!>  "empty", "unknown" and "error" results to the same return value (0),
!>  while ZSTD_getFrameContentSize() gives them separate return values.
!> @return : decompressed size of `src` frame content _if known and not empty_, 0 otherwise. */
!******************************************************************************!
        integer(c_long_long)function ZSTD_getDecompressedSize(src, srcSize) bind(C, name="ZSTD_getDecompressedSize")
            use iso_c_binding
            type(c_ptr), value :: src
            integer(c_size_t), value :: srcSize
        end function
!******************************************************************************!
!> `src` should point to the start of a ZSTD frame or skippable frame.
!> `srcSize` must be >= first frame size
!> @return : the compressed size of the first frame starting at `src`,
!>           suitable to pass as `srcSize` to `ZSTD_decompress` or similar,
!>           or an error code if input is invalid
!>@note
!> This method is called _find*() because it's not enough to read the header,
!> it may have to scan through the frame's content, to reach its end.
!>@endnote
!>
!>@note
!> This method also works with Skippable Frames. In which case,
!> it returns the size of the complete skippable frame,
!> which is always equal to its content size + 8 bytes for headers.
!>@endnote
!******************************************************************************!
        integer(c_long_long) function ZSTD_findFrameCompressedSize(src, srcSize) bind(C, name="ZSTD_findFrameCompressedSize")
            use iso_c_binding
            type(c_ptr), value :: src
            integer(c_size_t), value :: srcSize
        end function
!******************************************************************************!
!> maximum compressed size in worst case single-pass scenario.
!> When invoking `ZSTD_compress()`, or any other one-pass compression function,
!> it's recommended to provide @dstCapacity >= ZSTD_compressBound(srcSize)
!> as it eliminates one potential failure scenario,
!> aka not enough room in dst buffer to write the compressed frame.
!> Note : ZSTD_compressBound() itself can fail, if @srcSize >= ZSTD_MAX_INPUT_SIZE .
!>        In which case, ZSTD_compressBound() will return an error code
!>        which can be tested using ZSTD_isError().
!******************************************************************************!
        integer(c_size_t) function ZSTD_compressBound(srcSize) bind(C, name="ZSTD_compressBound")
            use iso_c_binding
            integer(c_size_t), value :: srcSize
        end function
!******************************************************************************!
!> Most ZSTD_* functions returning a size_t value can be tested for error,
!> using ZSTD_isError().
!> @return 1 if error, 0 otherwise
!******************************************************************************!
        logical(c_bool) function ZSTD_isError(code) bind(C, name="ZSTD_isError")
            use iso_c_binding
            integer(c_size_t), value :: code
        end function
!******************************************************************************!
!> Convert a result into an error code, which can be compared to error enum list
!******************************************************************************!
        integer(c_int) function ZSTD_getErrorCode(functionResult) bind(C, name="ZSTD_getErrorCode")
            use iso_c_binding
            integer(c_size_t), value :: functionResult
        end function
!******************************************************************************!
!> Provides readable string from a function result
!******************************************************************************!
        type(c_ptr) function ZSTD_getErrorName(code) bind(C, name="ZSTD_getErrorName")
            use iso_c_binding
            integer(c_size_t), value :: code
        end function
!******************************************************************************!
!> Minimum negative compression level allowed, requires v1.4.0+
!******************************************************************************!
        integer(c_int) function ZSTD_minCLevel() bind(C, name="ZSTD_minCLevel")
            use iso_c_binding
        end function
!******************************************************************************!
!> Maximum compression level available
!******************************************************************************!
        integer(c_int) function ZSTD_maxCLevel() bind(C, name="ZSTD_maxCLevel")
            use iso_c_binding
        end function
!******************************************************************************!
!> Default compression level, specified by `ZSTD_CLEVEL_DEFAULT`, requires v1.5.0+
!******************************************************************************!
        integer(c_int) function ZSTD_defaultCLevel() bind(C, name="ZSTD_defaultCLevel")
            use iso_c_binding
        end function
    end interface

    ! Explicit context
    interface
!******************************************************************************!
!>  When compressing many times,
!>  it is recommended to allocate a compression context just once,
!>  and reuse it for each successive compression operation.
!>  This will make the workload easier for system's memory.
!>  Note : re-using context is just a speed / resource optimization.
!>         It doesn't change the compression ratio, which remains identical.
!>  Note 2: For parallel execution in multi-threaded environments,
!>         use one different context per thread .
!******************************************************************************!
        type(ZSTD_CCtx) function ZSTD_createCCtx() bind(C, name="ZSTD_createCCtx")
            import
        end function
!******************************************************************************!
!******************************************************************************!
        integer(c_size_t) function ZSTD_freeCCtx(cctx) bind(C, name="ZSTD_freeCCtx")
            import
            type(ZSTD_CCtx), value :: cctx
        end function
!******************************************************************************!
!>  Same as ZSTD_compress(), using an explicit ZSTD_CCtx.
!>  Important : in order to mirror `ZSTD_compress()` behavior,
!>  this function compresses at the requested compression level,
!>  __ignoring any other advanced parameter__ .
!>  If any advanced parameter was set using the advanced API,
!>  they will all be reset. Only @compressionLevel remains.
!******************************************************************************!
        integer(c_size_t) function ZSTD_compressCCtx(cctx, dst, dstCapacity, src, srcSize, compressionLevel)&
            bind(C, name="ZSTD_compressCCtx")

            import
            type(ZSTD_CCtx), value :: cctx
            type(c_ptr), value :: dst
            integer(c_size_t), value :: dstCapacity
            type(c_ptr), value :: src
            integer(c_size_t), value :: srcSize
            integer(c_int), value :: compressionLevel
        end function
!******************************************************************************!
!>  When decompressing many times,
!>  it is recommended to allocate a context only once,
!>  and reuse it for each successive compression operation.
!>  This will make workload friendlier for system's memory.
!>  Use one context per thread for parallel execution.
!******************************************************************************!
        type(ZSTD_DCtx) function ZSTD_createDCtx() bind(C, name="ZSTD_createDCtx")
            import
        end function
!******************************************************************************!
!******************************************************************************!
        integer(c_size_t) function ZSTD_freeDCtx(cctx) bind(C, name="ZSTD_freeDCtx")
            import
            type(ZSTD_DCtx), value :: cctx
        end function
!******************************************************************************!
!> Same as [[ZSTD_decompress(function)]], requires an allocated `ZSTD_DCtx`.
!>
!> Compatible with sticky parameters.
!******************************************************************************!
        integer(c_size_t) function ZSTD_decompressDCtx(dctx, dst, dstCapacity, src, srcSize)&
            bind(C, name="ZSTD_decompressDCtx")

            import
            type(ZSTD_DCtx), value :: dctx
            type(c_ptr), value :: dst
            integer(c_size_t), value :: dstCapacity
            type(c_ptr), value :: src
            integer(c_size_t), value :: srcSize
        end function
    end interface

    type, public, bind(C) :: ZSTD_bounds
        integer(c_size_t) :: error
        integer(c_int) :: lowerBound
        integer(c_int) :: upperBound
    end type



    ! Streaming
    type, public, bind(C) :: ZSTD_inBuffer
        type(c_ptr) :: src = c_null_ptr !< start of input buffer
        integer(c_size_t) :: size = 0 !< size of input buffer
        integer(c_size_t) :: pos = 0 !< position where reading stopped. Will be updated. Necessarily 0 <= pos <= size
    end type

    type, public, bind(C) :: ZSTD_outBuffer
        type(c_ptr) :: dst = c_null_ptr !< start of output buffer
        integer(c_size_t) :: size = 0 !< size of output buffer
        integer(c_size_t) :: pos = 0 !< position where writing stopped. Will be updated. Necessarily 0 <= pos <= size
    end type

    enum, bind(C)
        !> Collect more data, encoder decides when to output compressed result, for optimal compression ratio
        enumerator :: ZSTD_e_continue = 0
        !> Flush any data provided so far, it creates (at least) one new block,
        !> that can be decoded immediately on reception; frame will continue:
        !> any future data can still reference previously compressed data,
        !> improving compression.
        !>@note
        !> Multithreaded compression will block to flush as much output as possible.
        !>@endnote
        enumerator :: ZSTD_e_flush    = 1
        !> Flush any remaining data _and_ close current frame. note that frame is only
        !> closed after compressed data is fully flushed (return value == 0). After that
        !> point, any additional data starts a new frame.
        !>@note
        !> Each frame is independent (does not reference any content from previous frame).
        !>@endnote
        !>@note
        !> multithreaded compression will block to flush as much output as possible.
        !>@endnot
        enumerator :: ZSTD_e_end      = 2
    end enum

    interface
!******************************************************************************!
!>  Behaves about the same as ZSTD_compressStream, with additional control on end directive.
!>  - Compression parameters are pushed into CCtx before starting compression, using ZSTD_CCtx_set*()
!>  - Compression parameters cannot be changed once compression is started (save a list of exceptions in multi-threading mode)
!>  - output->pos must be <= dstCapacity, input->pos must be <= srcSize
!>  - output->pos and input->pos will be updated. They are guaranteed to remain below their respective limit.
!>  - endOp must be a valid directive
!>  - When nbWorkers==0 (default), function is blocking : it completes its job before returning to caller.
!>  - When nbWorkers>=1, function is non-blocking : it copies a portion of input, distributes jobs to internal worker threads, flush to output whatever is available,
!>                                                  and then immediately returns, just indicating that there is some data remaining to be flushed.
!>                                                  The function nonetheless guarantees forward progress : it will return only after it reads or write at least 1+ byte.
!>  - Exception : if the first call requests a ZSTD_e_end directive and provides enough dstCapacity, the function delegates to ZSTD_compress2() which is always blocking.
!>  - @return provides a minimum amount of data remaining to be flushed from internal buffers
!>            or an error code, which can be tested using ZSTD_isError().
!>            if @return != 0, flush is not fully completed, there is still some data left within internal buffers.
!>            This is useful for ZSTD_e_flush, since in this case more flushes are necessary to empty all buffers.
!>            For ZSTD_e_end, @return == 0 when internal buffers are fully flushed and frame is completed.
!>  - after a ZSTD_e_end directive, if internal buffer is not fully flushed (@return != 0),
!>            only ZSTD_e_end or ZSTD_e_flush operations are allowed.
!>            Before starting a new compression job, or changing compression parameters,
!>            it is required to fully flush internal buffers.
!>  - note: if an operation ends with an error, it may leave @cctx in an undefined state.
!>          Therefore, it's UB to invoke ZSTD_compressStream2() of ZSTD_compressStream() on such a state.
!>          In order to be re-employed after an error, a state must be reset,
!>          which can be done explicitly (ZSTD_CCtx_reset()),
!>          or is sometimes implied by methods starting a new compression job (ZSTD_initCStream(), ZSTD_compressCCtx())
!******************************************************************************!
        integer(c_size_t) function ZSTD_compressStream2(cctx, output, input, endOp) bind(C, name="ZSTD_compressStream2")
            import
            type(ZSTD_CCtx), value :: cctx
            type(ZSTD_outBuffer), intent(inout) :: output
            type(ZSTD_inBuffer), intent(inout) :: input
            integer(c_int), value :: endOp
        end function
    end interface

    ! Legacy streaming APIs
    interface
        type(ZSTD_CStream) function ZSTD_createCStream() bind(C, name="ZSTD_createCStream")
            import
        end function
        integer(c_size_t) function ZSTD_freeCStream(zcs) bind(C, name="ZSTD_freeCStream")
            import
            type(ZSTD_CStream), value :: zcs
        end function
!******************************************************************************!
!> Equivalent to:
!>
!>     ZSTD_CCtx_reset(zcs, ZSTD_reset_session_only);
!>     ZSTD_CCtx_refCDict(zcs, NULL); // clear the dictionary (if any)
!>     ZSTD_CCtx_setParameter(zcs, ZSTD_c_compressionLevel, compressionLevel);
!>
!> Note that ZSTD_initCStream() clears any previously set dictionary. Use the new API
!> to compress with a dictionary.
!******************************************************************************!
        integer(c_size_t) function ZSTD_initCStream(zcs, compressionLevel) bind(C, name="ZSTD_initCStream")
            import
            type(ZSTD_CStream), value :: zcs
            integer(c_int), value :: compressionLevel
        end function
!******************************************************************************!
!> Alternative for `ZSTD_compressStream2(zcs, output, input, ZSTD_e_continue)`.
!>@note
!> The return value is different. ZSTD_compressStream() returns a hint for
!> the next read size (if non-zero and not an error). ZSTD_compressStream2()
!> returns the minimum nb of bytes left to flush (if non-zero and not an error).
!>@endnote
!******************************************************************************!
        integer(c_size_t) function ZSTD_compressStream(zcs, output, input) bind(C, name="ZSTD_compressStream")
            import
            type(ZSTD_CStream), value :: zcs
            type(ZSTD_outBuffer), intent(inout) :: output
            type(ZSTD_inBuffer), intent(inout) :: input
        end function
!******************************************************************************!
!> Equivalent to `ZSTD_compressStream2(zcs, output, &emptyInput, ZSTD_e_flush)`.
!******************************************************************************!
        integer(c_size_t) function ZSTD_flushStream(zcs, output) bind(C, name="ZSTD_flushStream")
            import
            type(ZSTD_CStream), value :: zcs
            type(ZSTD_outBuffer), intent(inout) :: output
        end function
!******************************************************************************!
!> Equivalent to `ZSTD_compressStream2(zcs, output, &emptyInput, ZSTD_e_end)`.
!******************************************************************************!
        integer(c_size_t) function ZSTD_endStream(zcs, output) bind(C, name="ZSTD_endStream")
            import
            type(ZSTD_CStream), value :: zcs
            type(ZSTD_outBuffer), intent(inout) :: output
        end function
    end interface

    interface
        type(ZSTD_DStream) function ZSTD_createDStream() bind(C, name="ZSTD_createDStream")
            import
        end function
        integer(c_size_t) function ZSTD_freeDStream(zds) bind(C, name="ZSTD_freeDStream")
            import
            type(ZSTD_DStream), value :: zds
        end function
!******************************************************************************!
!> Initialize/reset DStream state for new decompression operation.
!> Call before new decompression operation using same DStream.
!>
!>@note
!> This function is redundant with the advanced API and equivalent to:
!>     ZSTD_DCtx_reset(zds, ZSTD_reset_session_only);
!>     ZSTD_DCtx_refDDict(zds, NULL);
!>@endnote
!******************************************************************************!
        integer(c_size_t) function ZSTD_initDStream(zds) bind(C, name="ZSTD_initDStream")
            import
            type(ZSTD_DStream), value :: zds
        end function
!******************************************************************************!
!> Streaming decompression function.
!> Call repetitively to consume full input updating it as necessary.
!> Function will update both input and output `pos` fields exposing current state via these fields:
!> - `input.pos < input.size`, some input remaining and caller should provide remaining input
!>   on the next call.
!> - `output.pos < output.size`, decoder flushed internal output buffer.
!> - `output.pos == output.size`, unflushed data potentially present in the internal buffers,
!>   check ZSTD_decompressStream() @return value,
!>   if > 0, invoke it again to flush remaining data to output.
!>@note
!> With no additional input, amount of data flushed <= ZSTD_BLOCKSIZE_MAX.
!>@endnote
!>
!> @return : 0 when a frame is completely decoded and fully flushed,
!>           or an error code, which can be tested using ZSTD_isError(),
!>           or any other value > 0, which means there is some decoding or flushing to do to complete current frame.
!>
!>@note
!> When an operation returns with an error code, the @zds state may be left in undefined state.
!> It's UB to invoke `ZSTD_decompressStream()` on such a state. In order to re-use such a state,
!> it must be first reset, which can be done explicitly (`ZSTD_DCtx_reset()`),
!> or is implied for operations starting some new decompression job (`ZSTD_initDStream`,
!> `ZSTD_decompressDCtx()`, `ZSTD_decompress_usingDict()`)
!>@endnote
!******************************************************************************!
        integer(c_size_t) function ZSTD_decompressStream(zcs, output, input) bind(C, name="ZSTD_decompressStream")
            import
            type(ZSTD_DStream), value :: zcs
            type(ZSTD_outBuffer) :: output
            type(ZSTD_inBuffer) :: input
        end function
        integer(c_size_t) function ZSTD_DStreamInSize() bind(C, name="ZSTD_DStreamInSize")
            import
        end function
        integer(c_size_t) function ZSTD_DStreamOutSize() bind(C, name="ZSTD_DStreamOutSize")
            import
        end function
    end interface
end module
