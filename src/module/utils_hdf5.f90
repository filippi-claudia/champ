module hdf5_utils
    !> @details This module contains some utility functions for HDF5 data handeling.
    !> @author  Ravindra Shinde
    !> @email   r.l.shinde@utwente.nl
    !> @date    2022-11-03
    !> @version 1.0
    !> @note    This file is will be used to read and write HDF5 restart files.

    use hdf5
    use iso_c_binding,  only:   c_loc
    use contrl_file,    only:   errunit

    implicit none

    private
    public :: hdf5_file_create
    public :: hdf5_file_open
    public :: hdf5_file_close

    public :: hdf5_group_create
    public :: hdf5_group_open
    public :: hdf5_group_close

    public :: hdf5_write
    public :: hdf5_read

    public :: hid_t

    interface hdf5_write
        ! scalar variables
        module procedure hdf5_write_string
        module procedure hdf5_write_integer
        module procedure hdf5_write_real
        module procedure hdf5_write_double
        ! one dimensional arrays
        module procedure hdf5_write_array_integer_1d
        module procedure hdf5_write_array_real_1d
        module procedure hdf5_write_array_double_1d
        ! two dimensional arrays
        module procedure hdf5_write_array_integer_2d
        module procedure hdf5_write_array_real_2d
        module procedure hdf5_write_array_double_2d
        ! three dimensional arrays
        module procedure hdf5_write_array_integer_3d
        module procedure hdf5_write_array_real_3d
        module procedure hdf5_write_array_double_3d
        ! four dimensional arrays
        module procedure hdf5_write_array_integer_4d
        module procedure hdf5_write_array_real_4d
        module procedure hdf5_write_array_double_4d
    end interface hdf5_write

    interface hdf5_read
        ! scalar variables
        module procedure hdf5_read_string
        module procedure hdf5_read_integer
        module procedure hdf5_read_real
        module procedure hdf5_read_double
        ! ! one dimensional arrays
        module procedure hdf5_read_array_integer_1d
        module procedure hdf5_read_array_real_1d
        module procedure hdf5_read_array_double_1d
        ! ! two dimensional arrays
        module procedure hdf5_read_array_integer_2d
        module procedure hdf5_read_array_real_2d
        module procedure hdf5_read_array_double_2d
        ! ! three dimensional arrays
        module procedure hdf5_read_array_integer_3d
        module procedure hdf5_read_array_real_3d
        module procedure hdf5_read_array_double_3d
        ! ! four dimensional arrays
        module procedure hdf5_read_array_integer_4d
        module procedure hdf5_read_array_real_4d
        module procedure hdf5_read_array_double_4d
    end interface hdf5_read

    contains

    subroutine hdf5_file_create(filename, file_id)
        !> @brief Create a HDF5 file
        !> @param filename Name of the file to be created
        !> @param file_id  File id of the created file
        !> @author  Ravindra Shinde
        !> @email   r.l.shinde@utwente.nl
        use hdf5
        character(len=*), intent(in)        :: filename
        integer(hid_t), intent(out)         :: file_id
        integer                             :: ierr

        ! Initilize HDF5 Fortran
        call h5open_f(ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 library could not be initialized."
            stop
        end if

        ! create hdf5 file safely
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 file could not be created."
            stop
        end if

    end subroutine hdf5_file_create

    subroutine hdf5_file_open(filename, file_id)
        !> @brief Open a HDF5 file
        !> @param filename Name of the file to be opened
        !> @param file_id  File id of the opened file
        !> @author  Ravindra Shinde
        !> @email   r.l.shinde@utwente.nl

        character(len=*), intent(in)        :: filename
        integer(hid_t), intent(out)         :: file_id
        integer                             :: ierr

        ! Initilize HDF5 Fortran
        call h5open_f(ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 library could not be initialized."
            stop
        end if

        ! open file for reading and writing
        call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 file could not be opened."
            stop
        end if

    end subroutine hdf5_file_open

    subroutine hdf5_file_close(file_id)
        !> @brief Close a HDF5 file
        !> @param file_id  File id of the file to be closed
        !> @author  Ravindra Shinde
        !> @email   r.l.shinde@utwente.nl
        integer(hid_t), intent(inout)          :: file_id
        integer                             :: ierr

        ! close file
        call h5fclose_f(file_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 file could not be closed."
            stop
        end if

        ! close HDF5 Fortran
        call h5close_f(ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 library could not be closed."
            stop
        end if

    end subroutine hdf5_file_close

    subroutine hdf5_group_create(file_id, group_name, group_id)
        !> @brief Create a HDF5 group
        !> @param file_id    File id of the file in which the group is created
        !> @param group_name Name of the group to be created
        !> @param group_id   Group id of the created group
        !> @author  Ravindra Shinde
        !> @email   r.l.shinde@utwente.nl

        integer(hid_t), intent(in)             :: file_id
        character(len=*), intent(in)        :: group_name
        integer(hid_t), intent(out)            :: group_id
        integer                             :: ierr

        ! create group
        call h5gcreate_f(file_id, group_name, group_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 group could not be created."
            stop
        end if

    end subroutine hdf5_group_create

    subroutine hdf5_group_open(file_id, group_name, group_id)
        !> @brief Open a HDF5 group
        !> @param file_id    File id of the file in which the group is opened
        !> @param group_name Name of the group to be opened
        !> @param group_id   Group id of the opened group
        !> @author  Ravindra Shinde
        !> @email   r.l.shinde@utwente.nl

        integer(hid_t), intent(in)             :: file_id
        character(len=*), intent(in)        :: group_name
        integer(hid_t), intent(out)            :: group_id
        integer                             :: ierr

        ! open group
        call h5gopen_f(file_id, group_name, group_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 group could not be opened."
            stop
        end if

    end subroutine hdf5_group_open

    subroutine hdf5_group_close(group_id)
        !> @brief Close a HDF5 group
        !> @param group_id   Group id of the group to be closed
        !> @author  Ravindra Shinde
        !> @email   r.l.shinde@utwente.nl

        integer(hid_t), intent(inout)          :: group_id
        integer                             :: ierr

        ! close group
        call h5gclose_f(group_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 group could not be closed."
            stop
        end if

    end subroutine hdf5_group_close

    ! The following subroutine writes integer number to a HDF5 file
    subroutine hdf5_write_integer(file_id, group_id, dataset_name, data)
        !> @brief Write a integer number to a HDF5 file
        !> @param file_id       File id of the file in which the data is written
        !> @param group_id      Group id in which the data is written
        !> @param dataset_name  Name of the dataset in which the data is written
        !> @param data          Data to be written
        !> @author  Ravindra Shinde
        !> @email   r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        integer, intent(in)                     :: data
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer                                 :: ierr


        ! create dataspace
        call h5screate_f(H5S_SCALAR_F, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be created."
            stop
        end if

        ! create dataset
        call h5dcreate_f(group_id, dataset_name, H5T_NATIVE_INTEGER, dataspace_id, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be created."
            stop
        end if

        ! write data
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, data, [0_HSIZE_T,0_HSIZE_T], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be written."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_write_integer


    ! The following subroutine writes string to a HDF5 file
    subroutine hdf5_write_string(file_id, group_id, dataset_name, data)
        !> @brief Write a string to a HDF5 file
        !> @param file_id       File id of the file in which the data is written
        !> @param group_id      Group id in which the data is written
        !> @param dataset_name  Name of the dataset in which the data is written
        !> @param data          Data to be written
        !> @author  Ravindra Shinde
        !> @email   r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        character(len=*), intent(in)            :: data
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(hid_t)                          :: type_id
        integer                                 :: ierr


        ! get length of data string
        call h5tcopy_f(H5T_FORTRAN_S1, type_id, ierr)
        call h5tset_size_f(type_id, len(data, SIZE_T), ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 could not obtain size of a string."
            stop
        end if


        ! create dataspace
        call h5screate_f(H5S_SCALAR_F, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be created."
            stop
        end if

        ! create dataset
        call h5dcreate_f(group_id, dataset_name, type_id, dataspace_id, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be created."
            stop
        end if

        ! write data
        call h5dwrite_f(dataset_id, type_id, data, [0_HSIZE_T], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be written."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*), "Error: HDF5 dataspace could not be closed."
            stop
        end if

        ! Release fortran string.
        call h5tclose_f(type_id, ierr)

    end subroutine hdf5_write_string


   ! The following subroutine writes real number to a HDF5 file
    subroutine hdf5_write_real(file_id, group_id, dataset_name, data)
        !> @brief Write a real number to a HDF5 file
        !> @param file_id       File id of the file in which the data is written
        !> @param group_id      Group id in which the data is written
        !> @param dataset_name  Name of the dataset in which the data is written
        !> @param data          Data to be written
        !> @author  Ravindra Shinde
        !> @email   r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        real, intent(in)                        :: data
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer                                 :: ierr


        ! create dataspace
        call h5screate_f(H5S_SCALAR_F, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be created."
            stop
        end if

        ! create dataset
        call h5dcreate_f(group_id, dataset_name, H5T_NATIVE_REAL, dataspace_id, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be created."
            stop
        end if

        ! write data
        call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, [0_HSIZE_T,0_HSIZE_T], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be written."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_write_real

    ! The following subroutine writes a double precision number to a HDF5 file
    subroutine hdf5_write_double(file_id, group_id, dataset_name, data)
        !> @brief Write a double precision number to a HDF5 file
        !> @param file_id       File id of the file in which the data is written
        !> @param group_id      Group id in which the data is written
        !> @param dataset_name  Name of the dataset in which the data is written
        !> @param data          Data to be written
        !> @author  Ravindra Shinde
        !> @email   r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        double precision, intent(in)            :: data
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer                                 :: ierr

        ! create dataspace
        call h5screate_f(H5S_SCALAR_F, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be created."
            stop
        end if

        ! create dataset
        call h5dcreate_f(group_id, dataset_name, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be created."
            stop
        end if

        ! write data
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, [0_HSIZE_T,0_HSIZE_T], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be written."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_write_double

    ! The following subroutine writes real array to a HDF5 file
    subroutine hdf5_write_array_real_1d(file_id, group_id, dataset_name, data)
        !> @brief Write a real 1d array to a HDF5 file
        !> @param file_id       File id of the file in which the data is written
        !> @param group_id      Group id in which the data is written
        !> @param dataset_name  Name of the dataset in which the data is written
        !> @param data          Data to be written
        !> @author  Ravindra Shinde
        !> @email   r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        real, intent(in)                        :: data(:)
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_size
        integer(HSIZE_T)                        :: data_dims = 1
        integer                                 :: ierr

        data_size = size(data)

        ! create dataspace
        call h5screate_f(H5S_SIMPLE_F, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be created."
            stop
        end if

        ! set dimensions
        call h5sset_extent_simple_f(dataspace_id, 1, [data_size], [data_size], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be set."
            stop
        end if

        ! create dataset
        call h5dcreate_f(group_id, dataset_name, H5T_NATIVE_REAL, dataspace_id, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be created."
            stop
        end if

        ! write data
        call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, [0_HSIZE_T,0_HSIZE_T], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be written."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_write_array_real_1d


    ! The following subroutine writes a real 2D array to a HDF5 file
    subroutine hdf5_write_array_real_2d(file_id, group_id, dataset_name, data)
        !> @brief Write a real 2D array to a HDF5 file
        !> @param file_id       File id of the file in which the data is written
        !> @param group_id      Group id in which the data is written
        !> @param dataset_name  Name of the dataset in which the data is written
        !> @param data          Data to be written
        !> @author  Ravindra Shinde
        !> @email   r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        real, intent(in)                        :: data(:,:)
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_shape(2)
        integer(HSIZE_T)                        :: data_dims = 2
        integer                                 :: ierr

        data_shape = shape(data)

        ! create dataspace
        call h5screate_f(H5S_SIMPLE_F, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be created."
            stop
        end if

        ! set dimensions
        call h5sset_extent_simple_f(dataspace_id, 2, [data_shape(1), data_shape(2)], [data_shape(1), data_shape(2)], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be set."
            stop
        end if

        ! create dataset
        call h5dcreate_f(group_id, dataset_name, H5T_NATIVE_REAL, dataspace_id, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be created."
            stop
        end if

        ! write data
        call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, [0_HSIZE_T,0_HSIZE_T], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be written."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_write_array_real_2d


    ! The following subroutine writes a real 3D array to a HDF5 file
    subroutine hdf5_write_array_real_3d(file_id, group_id, dataset_name, data)
        !> @brief Write a real 3D array to a HDF5 file
        !> @param file_id       File id of the file in which the data is written
        !> @param group_id      Group id in which the data is written
        !> @param dataset_name  Name of the dataset in which the data is written
        !> @param data          Data to be written
        !> @author  Ravindra Shinde
        !> @email   r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        real, intent(in)                        :: data(:,:,:)
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_shape(3)
        integer(HSIZE_T)                        :: data_dims = 3
        integer                                 :: ierr

        data_shape = shape(data)

        ! create dataspace
        call h5screate_f(H5S_SIMPLE_F, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be created."
            stop
        end if

        ! set dimensions
        call h5sset_extent_simple_f(dataspace_id, 3, [data_shape(1), data_shape(2), data_shape(3)], [data_shape(1), data_shape(2), data_shape(3)], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be set."
            stop
        end if

        ! create dataset
        call h5dcreate_f(group_id, dataset_name, H5T_NATIVE_REAL, dataspace_id, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be created."
            stop
        end if

        ! write data
        call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, [0_HSIZE_T,0_HSIZE_T,0_HSIZE_T], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be written."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_write_array_real_3d


    ! The following subroutine writes a real 4D array to a HDF5 file
    subroutine hdf5_write_array_real_4d(file_id, group_id, dataset_name, data)
        !> @brief Write a real 4D array to a HDF5 file
        !> @param file_id       File id of the file in which the data is written
        !> @param group_id      Group id in which the data is written
        !> @param dataset_name  Name of the dataset in which the data is written
        !> @param data          Data to be written
        !> @author  Ravindra Shinde
        !> @email   r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        real, intent(in)                        :: data(:,:,:,:)
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_shape(4)
        integer(HSIZE_T)                        :: data_dims = 4
        integer                                 :: ierr

        data_shape = shape(data)

        ! create dataspace
        call h5screate_f(H5S_SIMPLE_F, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be created."
            stop
        end if

        ! set dimensions
        call h5sset_extent_simple_f(dataspace_id, 4, [data_shape(1), data_shape(2), data_shape(3), data_shape(4)], [data_shape(1), data_shape(2), data_shape(3), data_shape(4)], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be set."
            stop
        end if

        ! create dataset
        call h5dcreate_f(group_id, dataset_name, H5T_NATIVE_REAL, dataspace_id, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be created."
            stop
        end if

        ! write data
        call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, [0_HSIZE_T,0_HSIZE_T,0_HSIZE_T,0_HSIZE_T], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be written."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_write_array_real_4d

    ! The following subroutine writes a double precision 1D array to a HDF5 file
    subroutine hdf5_write_array_double_1d(file_id, group_id, dataset_name, data)
        !> @brief Write a double precision 1D array to a HDF5 file
        !> @param file_id       File id of the file in which the data is written
        !> @param group_id      Group id in which the data is written
        !> @param dataset_name  Name of the dataset in which the data is written
        !> @param data          Data to be written
        !> @author  Ravindra Shinde
        !> @email   r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        double precision, intent(in)            :: data(:)
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_shape(1)
        integer(HSIZE_T)                        :: data_dims = 1
        integer                                 :: ierr

        data_shape = shape(data)

        ! create dataspace
        call h5screate_f(H5S_SIMPLE_F, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be created."
            stop
        end if

        ! set dimensions
        call h5sset_extent_simple_f(dataspace_id, 1, [data_shape(1)], [data_shape(1)], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be set."
            stop
        end if

        ! create dataset
        call h5dcreate_f(group_id, dataset_name, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be created."
            stop
        end if

        ! write data
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, [0_HSIZE_T], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be written."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_write_array_double_1d


    ! The following subroutine writes a double precision 2D array to a HDF5 file
    subroutine hdf5_write_array_double_2d(file_id, group_id, dataset_name, data)
        !> @brief Write a double precision 2D array to a HDF5 file
        !> @param file_id       File id of the file in which the data is written
        !> @param group_id      Group id in which the data is written
        !> @param dataset_name  Name of the dataset in which the data is written
        !> @param data          Data to be written
        !> @author  Ravindra Shinde
        !> @email   r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        double precision, intent(in)            :: data(:,:)
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_shape(2)
        integer(HSIZE_T)                        :: data_dims = 2
        integer                                 :: ierr

        data_shape = shape(data)

        ! create dataspace
        call h5screate_f(H5S_SIMPLE_F, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be created."
            stop
        end if

        ! set dimensions
        call h5sset_extent_simple_f(dataspace_id, 2, [data_shape(1), data_shape(2)], [data_shape(1), data_shape(2)], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be set."
            stop
        end if

        ! create dataset
        call h5dcreate_f(group_id, dataset_name, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be created."
            stop
        end if

        ! write data
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, [0_HSIZE_T,0_HSIZE_T], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be written."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_write_array_double_2d

    ! The following subroutine writes a double precision 3D array to a HDF5 file
    subroutine hdf5_write_array_double_3d(file_id, group_id, dataset_name, data)
        !> @brief Write a double precision 3D array to a HDF5 file
        !> @param file_id       File id of the file in which the data is written
        !> @param group_id      Group id in which the data is written
        !> @param dataset_name  Name of the dataset in which the data is written
        !> @param data          Data to be written
        !> @author  Ravindra Shinde
        !> @email   r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        double precision, intent(in)            :: data(:,:,:)
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_shape(3)
        integer(HSIZE_T)                        :: data_dims = 3
        integer                                 :: ierr

        data_shape = shape(data)

        ! create dataspace
        call h5screate_f(H5S_SIMPLE_F, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be created."
            stop
        end if

        ! set dimensions
        call h5sset_extent_simple_f(dataspace_id, 3, [data_shape(1), data_shape(2), data_shape(3)], [data_shape(1), data_shape(2), data_shape(3)], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be set."
            stop
        end if

        ! create dataset
        call h5dcreate_f(group_id, dataset_name, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be created."
            stop
        end if

        ! write data
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, [0_HSIZE_T,0_HSIZE_T,0_HSIZE_T], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be written."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_write_array_double_3d

    ! The following subroutine writes a double precision 4D array to a HDF5 file
    subroutine hdf5_write_array_double_4d(file_id, group_id, dataset_name, data)
        !> @brief Write a double precision 4D array to a HDF5 file
        !> @param file_id       File id of the file in which the data is written
        !> @param group_id      Group id in which the data is written
        !> @param dataset_name  Name of the dataset in which the data is written
        !> @param data          Data to be written
        !> @author  Ravindra Shinde
        !> @email   r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        double precision, intent(in)            :: data(:,:,:,:)
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_shape(4)
        integer(HSIZE_T)                        :: data_dims = 4
        integer                                 :: ierr

        data_shape = shape(data)

        ! create dataspace
        call h5screate_f(H5S_SIMPLE_F, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be created."
            stop
        end if

        ! set dimensions
        call h5sset_extent_simple_f(dataspace_id, 4, [data_shape(1), data_shape(2), data_shape(3), data_shape(4)], [data_shape(1), data_shape(2), data_shape(3), data_shape(4)], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be set."
            stop
        end if

        ! create dataset
        call h5dcreate_f(group_id, dataset_name, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be created."
            stop
        end if

        ! write data
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, [0_HSIZE_T,0_HSIZE_T,0_HSIZE_T,0_HSIZE_T], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be written."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_write_array_double_4d

    ! The following subroutine writes an integer array to a HDF5 file
    subroutine hdf5_write_array_integer_1d(file_id, group_id, dataset_name, data)
        !> @brief Write an integer array to a HDF5 file
        !> @param file_id       File id of the file in which the data is written
        !> @param group_id      Group id in which the data is written
        !> @param dataset_name  Name of the dataset in which the data is written
        !> @param data          Data to be written
        !> @author  Ravindra Shinde
        !> @email   r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        integer, intent(in)                     :: data(:)
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_size
        integer(HSIZE_T)                        :: data_dims = 1
        integer                                 :: ierr

        data_size = size(data)

        ! create dataspace
        call h5screate_f(H5S_SIMPLE_F, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be created."
            stop
        end if

        ! set dimensions
        call h5sset_extent_simple_f(dataspace_id, 1, [data_size], [data_size], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be set."
            stop
        end if

        ! create dataset
        call h5dcreate_f(group_id, dataset_name, H5T_NATIVE_INTEGER, dataspace_id, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be created."
            stop
        end if

        ! write data
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, data, [0_HSIZE_T,0_HSIZE_T], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be written."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_write_array_integer_1d


    ! The following subroutine writes an integer 2D array to a HDF5 file
    subroutine hdf5_write_array_integer_2d(file_id, group_id, dataset_name, data)
        !> @brief Write an integer 2D array to a HDF5 file
        !> @param file_id       File id of the file in which the data is written
        !> @param group_id      Group id in which the data is written
        !> @param dataset_name  Name of the dataset in which the data is written
        !> @param data          Data to be written
        !> @author  Ravindra Shinde
        !> @email   r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        integer, intent(in)                     :: data(:,:)
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_shape(2)
        integer(HSIZE_T)                        :: data_dims = 2
        integer                                 :: ierr

        data_shape = shape(data)

        ! create dataspace
        call h5screate_f(H5S_SIMPLE_F, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be created."
            stop
        end if

        ! set dimensions
        call h5sset_extent_simple_f(dataspace_id, 2, [data_shape(1), data_shape(2)], [data_shape(1), data_shape(2)], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be set."
            stop
        end if

        ! create dataset
        call h5dcreate_f(group_id, dataset_name, H5T_NATIVE_INTEGER, dataspace_id, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be created."
            stop
        end if

        ! write data
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, data, [0_HSIZE_T,0_HSIZE_T], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be written."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_write_array_integer_2d

    ! The following subroutine writes an integer 3D array to a HDF5 file
    subroutine hdf5_write_array_integer_3d(file_id, group_id, dataset_name, data)
        !> @brief Write an integer 3D array to a HDF5 file
        !> @param file_id       File id of the file in which the data is written
        !> @param group_id      Group id in which the data is written
        !> @param dataset_name  Name of the dataset in which the data is written
        !> @param data          Data to be written
        !> @author  Ravindra Shinde
        !> @email   r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        integer, intent(in)                     :: data(:,:,:)
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_shape(3)
        integer(HSIZE_T)                        :: data_dims = 3
        integer                                 :: ierr

        data_shape = shape(data)

        ! create dataspace
        call h5screate_f(H5S_SIMPLE_F, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be created."
            stop
        end if

        ! set dimensions
        call h5sset_extent_simple_f(dataspace_id, 3, [data_shape(1), data_shape(2), data_shape(3)], [data_shape(1), data_shape(2), data_shape(3)], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be set."
            stop
        end if

        ! create dataset
        call h5dcreate_f(group_id, dataset_name, H5T_NATIVE_INTEGER, dataspace_id, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be created."
            stop
        end if

        ! write data
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, data, [0_HSIZE_T,0_HSIZE_T,0_HSIZE_T], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be written."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_write_array_integer_3d

    ! The following subroutine writes an integer 4D array to a HDF5 file
    subroutine hdf5_write_array_integer_4d(file_id, group_id, dataset_name, data)
        !> @brief Write an integer 4D array to a HDF5 file
        !> @param file_id       File id of the file in which the data is written
        !> @param group_id      Group id in which the data is written
        !> @param dataset_name  Name of the dataset in which the data is written
        !> @param data          Data to be written
        !> @author  Ravindra Shinde
        !> @email   r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        integer, intent(in)                     :: data(:,:,:,:)
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_shape(4)
        integer(HSIZE_T)                        :: data_dims = 4
        integer                                 :: ierr

        data_shape = shape(data)

        ! create dataspace
        call h5screate_f(H5S_SIMPLE_F, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be created."
            stop
        end if

        ! set dimensions
        call h5sset_extent_simple_f(dataspace_id, 4, [data_shape(1), data_shape(2), data_shape(3), data_shape(4)], [data_shape(1), data_shape(2), data_shape(3), data_shape(4)], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be set."
            stop
        end if

        ! create dataset
        call h5dcreate_f(group_id, dataset_name, H5T_NATIVE_INTEGER, dataspace_id, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be created."
            stop
        end if

        ! write data
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, data, [0_HSIZE_T,0_HSIZE_T,0_HSIZE_T,0_HSIZE_T], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be written."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_write_array_integer_4d


    ! The Following ste of subroutines are used for reading the data from a HDF5 file

    ! The following subroutine reads an integer from a HDF5 file
    subroutine hdf5_read_integer(file_id, group_id, dataset_name, data)
        !> @brief Read an integer from a HDF5 file
        !> @param file_id       File id of the file from which the data is read
        !> @param group_id      Group id from which the data is read
        !> @param dataset_name  Name of the dataset from which the data is read
        !> @param data          Data to be read
        !> @author  Ravindra Shinde
        !> @email r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        integer, intent(out)                    :: data
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_dims
        integer(HSIZE_T)                        :: data_shape(1)
        integer                                 :: ierr

        ! open dataset
        call h5dopen_f(group_id, dataset_name, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be opened."
            stop
        end if

        ! get dataspace
        call h5dget_space_f(dataset_id, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be obtained."
            stop
        end if

        ! read data
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, data, [0_HSIZE_T], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be read."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_read_integer

    ! The following subroutine reads a real from a HDF5 file
    subroutine hdf5_read_real(file_id, group_id, dataset_name, data)
        !> @brief Read a real from a HDF5 file
        !> @param file_id       File id of the file from which the data is read
        !> @param group_id      Group id from which the data is read
        !> @param dataset_name  Name of the dataset from which the data is read
        !> @param data          Data to be read
        !> @author  Ravindra Shinde
        !> @email r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        real, intent(out)                       :: data
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_dims
        integer(HSIZE_T)                        :: data_shape(1)
        integer                                 :: ierr

        ! open dataset
        call h5dopen_f(group_id, dataset_name, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be opened."
            stop
        end if

        ! get dataspace
        call h5dget_space_f(dataset_id, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be obtained."
            stop
        end if

        ! read data
        call h5dread_f(dataset_id, H5T_NATIVE_REAL, data, [0_HSIZE_T], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be read."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_read_real

    ! The following subroutine reads a double precision real from a HDF5 file
    subroutine hdf5_read_double(file_id, group_id, dataset_name, data)
        !> @brief Read a double precision real from a HDF5 file
        !> @param file_id       File id of the file from which the data is read
        !> @param group_id      Group id from which the data is read
        !> @param dataset_name  Name of the dataset from which the data is read
        !> @param data          Data to be read
        !> @author  Ravindra Shinde
        !> @email r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        real(kind=8), intent(out)               :: data
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_dims
        integer(HSIZE_T)                        :: data_shape(1)
        integer                                 :: ierr

        ! open dataset
        call h5dopen_f(group_id, dataset_name, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be opened."
            stop
        end if

        ! get dataspace
        call h5dget_space_f(dataset_id, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be obtained."
            stop
        end if

        ! read data
        call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, [0_HSIZE_T], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be read."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_read_double


    ! The following subroutine reads a string from a HDF5 file
    subroutine hdf5_read_string(file_id, group_id, dataset_name, data)
        !> @brief Read a string from a HDF5 file
        !> @param file_id       File id of the file from which the data is read
        !> @param group_id      Group id from which the data is read
        !> @param dataset_name  Name of the dataset from which the data is read
        !> @param data          Data to be read
        !> @author  Ravindra Shinde
        !> @email r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        character(len=*), intent(out)           :: data
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_dims
        integer(HSIZE_T)                        :: data_shape(1)
        integer                                 :: ierr

        ! open dataset
        call h5dopen_f(group_id, dataset_name, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be opened."
            stop
        end if

        ! get dataspace
        call h5dget_space_f(dataset_id, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be obtained."
            stop
        end if

        ! read data
        call h5dread_f(dataset_id, H5T_C_S1, data, [0_HSIZE_T], ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be read."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_read_string

    ! The following subroutine reads a 1D array of integers from a HDF5 file
    subroutine hdf5_read_array_integer_1d(file_id, group_id, dataset_name, data)
        !> @brief Read a 1D array of integers from a HDF5 file
        !> @param file_id       File id of the file from which the data is read
        !> @param group_id      Group id from which the data is read
        !> @param dataset_name  Name of the dataset from which the data is read
        !> @param data          Data to be read
        !> @author  Ravindra Shinde
        !> @email r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        integer, intent(out)                    :: data(:)
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_dims
        integer(HSIZE_T)                        :: data_shape(1)
        integer                                 :: ierr

        data_shape = shape(data)

        ! open dataset
        call h5dopen_f(group_id, dataset_name, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be opened."
            stop
        end if

        ! get dataspace
        call h5dget_space_f(dataset_id, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be obtained."
            stop
        end if

        ! read data
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, data, data_shape, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be read."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_read_array_integer_1d

    ! The following subroutine reads a 2D array of integers from a HDF5 file
    subroutine hdf5_read_array_integer_2d(file_id, group_id, dataset_name, data)
        !> @brief Read a 2D array of integers from a HDF5 file
        !> @param file_id       File id of the file from which the data is read
        !> @param group_id      Group id from which the data is read
        !> @param dataset_name  Name of the dataset from which the data is read
        !> @param data          Data to be read
        !> @author  Ravindra Shinde
        !> @email r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        integer, intent(out)                    :: data(:,:)
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_dims
        integer(HSIZE_T)                        :: data_shape(2)
        integer                                 :: ierr

        data_shape = shape(data)

        ! open dataset
        call h5dopen_f(group_id, dataset_name, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be opened."
            stop
        end if

        ! get dataspace
        call h5dget_space_f(dataset_id, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be obtained."
            stop
        end if

        ! read data
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, data, data_shape, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be read."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_read_array_integer_2d

    ! The following subroutine reads a 3D array of integers from a HDF5 file
    subroutine hdf5_read_array_integer_3d(file_id, group_id, dataset_name, data)
        !> @brief Read a 3D array of integers from a HDF5 file
        !> @param file_id       File id of the file from which the data is read
        !> @param group_id      Group id from which the data is read
        !> @param dataset_name  Name of the dataset from which the data is read
        !> @param data          Data to be read
        !> @author  Ravindra Shinde
        !> @email r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        integer, intent(out)                    :: data(:,:,:)
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_dims
        integer(HSIZE_T)                        :: data_shape(3)
        integer                                 :: ierr

        data_shape = shape(data)

        ! open dataset
        call h5dopen_f(group_id, dataset_name, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be opened."
            stop
        end if

        ! get dataspace
        call h5dget_space_f(dataset_id, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be obtained."
            stop
        end if

        ! read data
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, data, data_shape, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be read."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_read_array_integer_3d

    ! The following subroutine reads a 4D array of integers from a HDF5 file
    subroutine hdf5_read_array_integer_4d(file_id, group_id, dataset_name, data)
        !> @brief Read a 4D array of integers from a HDF5 file
        !> @param file_id       File id of the file from which the data is read
        !> @param group_id      Group id from which the data is read
        !> @param dataset_name  Name of the dataset from which the data is read
        !> @param data          Data to be read
        !> @author  Ravindra Shinde
        !> @email r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        integer, intent(out)                    :: data(:,:,:,:)
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_dims
        integer(HSIZE_T)                        :: data_shape(4)
        integer                                 :: ierr

        data_shape = shape(data)

        ! open dataset
        call h5dopen_f(group_id, dataset_name, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be opened."
            stop
        end if

        ! get dataspace
        call h5dget_space_f(dataset_id, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be obtained."
            stop
        end if

        ! read data
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, data, data_shape, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be read."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_read_array_integer_4d

    ! The following subroutine reads a 1D array of reals from a HDF5 file
    subroutine hdf5_read_array_real_1d(file_id, group_id, dataset_name, data)
        !> @brief Read a 1D array of reals from a HDF5 file
        !> @param file_id       File id of the file from which the data is read
        !> @param group_id      Group id from which the data is read
        !> @param dataset_name  Name of the dataset from which the data is read
        !> @param data          Data to be read
        !> @author  Ravindra Shinde
        !> @email r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        real, intent(out)                       :: data(:)
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_dims
        integer(HSIZE_T)                        :: data_shape(1)
        integer                                 :: ierr

        data_shape = shape(data)

        ! open dataset
        call h5dopen_f(group_id, dataset_name, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be opened."
            stop
        end if

        ! get dataspace
        call h5dget_space_f(dataset_id, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be obtained."
            stop
        end if

        ! read data
        call h5dread_f(dataset_id, H5T_NATIVE_REAL, data, data_shape, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be read."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_read_array_real_1d

    ! The following subroutine reads a 2D array of reals from a HDF5 file
    subroutine hdf5_read_array_real_2d(file_id, group_id, dataset_name, data)
        !> @brief Read a 2D array of reals from a HDF5 file
        !> @param file_id       File id of the file from which the data is read
        !> @param group_id      Group id from which the data is read
        !> @param dataset_name  Name of the dataset from which the data is read
        !> @param data          Data to be read
        !> @author  Ravindra Shinde
        !> @email r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        real, intent(out)                       :: data(:,:)
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_dims
        integer(HSIZE_T)                        :: data_shape(2)
        integer                                 :: ierr

        data_shape = shape(data)

        ! open dataset
        call h5dopen_f(group_id, dataset_name, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be opened."
            stop
        end if

        ! get dataspace
        call h5dget_space_f(dataset_id, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be obtained."
            stop
        end if

        ! read data
        call h5dread_f(dataset_id, H5T_NATIVE_REAL, data, data_shape, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be read."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_read_array_real_2d

    ! The following subroutine reads a 3D array of reals from a HDF5 file
    subroutine hdf5_read_array_real_3d(file_id, group_id, dataset_name, data)
        !> @brief Read a 3D array of reals from a HDF5 file
        !> @param file_id       File id of the file from which the data is read
        !> @param group_id      Group id from which the data is read
        !> @param dataset_name  Name of the dataset from which the data is read
        !> @param data          Data to be read
        !> @author  Ravindra Shinde
        !> @email r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        real, intent(out)                       :: data(:,:,:)
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_dims
        integer(HSIZE_T)                        :: data_shape(3)
        integer                                 :: ierr

        data_shape = shape(data)

        ! open dataset
        call h5dopen_f(group_id, dataset_name, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be opened."
            stop
        end if

        ! get dataspace
        call h5dget_space_f(dataset_id, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be obtained."
            stop
        end if

        ! read data
        call h5dread_f(dataset_id, H5T_NATIVE_REAL, data, data_shape, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be read."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_read_array_real_3d

    ! The following subroutine reads a 4D array of reals from a HDF5 file
    subroutine hdf5_read_array_real_4d(file_id, group_id, dataset_name, data)
        !> @brief Read a 4D array of reals from a HDF5 file
        !> @param file_id       File id of the file from which the data is read
        !> @param group_id      Group id from which the data is read
        !> @param dataset_name  Name of the dataset from which the data is read
        !> @param data          Data to be read
        !> @author  Ravindra Shinde
        !> @email r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        real, intent(out)                       :: data(:,:,:,:)
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_dims
        integer(HSIZE_T)                        :: data_shape(4)
        integer                                 :: ierr

        data_shape = shape(data)

        ! open dataset
        call h5dopen_f(group_id, dataset_name, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be opened."
            stop
        end if

        ! get dataspace
        call h5dget_space_f(dataset_id, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be obtained."
            stop
        end if

        ! read data
        call h5dread_f(dataset_id, H5T_NATIVE_REAL, data, data_shape, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be read."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_read_array_real_4d

    ! The following subroutine reads a 1D array of double precision reals from a HDF5 file
    subroutine hdf5_read_array_double_1d(file_id, group_id, dataset_name, data)
        !> @brief Read a 1D array of double precision reals from a HDF5 file
        !> @param file_id       File id of the file from which the data is read
        !> @param group_id      Group id from which the data is read
        !> @param dataset_name  Name of the dataset from which the data is read
        !> @param data          Data to be read
        !> @author  Ravindra Shinde
        !> @email r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        double precision, intent(out)           :: data(:)
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_dims
        integer(HSIZE_T)                        :: data_shape(1)
        integer                                 :: ierr

        data_shape = shape(data)

        ! open dataset
        call h5dopen_f(group_id, dataset_name, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be opened."
            stop
        end if

        ! get dataspace
        call h5dget_space_f(dataset_id, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be obtained."
            stop
        end if

        ! read data
        call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, data_shape, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be read."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_read_array_double_1d

    ! The following subroutine reads a 2D array of double precision reals from a HDF5 file
    subroutine hdf5_read_array_double_2d(file_id, group_id, dataset_name, data)
        !> @brief Read a 2D array of double precision reals from a HDF5 file
        !> @param file_id       File id of the file from which the data is read
        !> @param group_id      Group id from which the data is read
        !> @param dataset_name  Name of the dataset from which the data is read
        !> @param data          Data to be read
        !> @author  Ravindra Shinde
        !> @email r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        double precision, intent(out)           :: data(:,:)
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_dims
        integer(HSIZE_T)                        :: data_shape(2)
        integer                                 :: ierr

        data_shape = shape(data)

        ! open dataset
        call h5dopen_f(group_id, dataset_name, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be opened."
            stop
        end if

        ! get dataspace
        call h5dget_space_f(dataset_id, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be obtained."
            stop
        end if

        ! read data
        call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, data_shape, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be read."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_read_array_double_2d

    ! The following subroutine reads a 3D array of double precision reals from a HDF5 file
    subroutine hdf5_read_array_double_3d(file_id, group_id, dataset_name, data)
        !> @brief Read a 3D array of double precision reals from a HDF5 file
        !> @param file_id       File id of the file from which the data is read
        !> @param group_id      Group id from which the data is read
        !> @param dataset_name  Name of the dataset from which the data is read
        !> @param data          Data to be read
        !> @author  Ravindra Shinde
        !> @email r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        double precision, intent(out)           :: data(:,:,:)
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_dims
        integer(HSIZE_T)                        :: data_shape(3)
        integer                                 :: ierr

        data_shape = shape(data)

        ! open dataset
        call h5dopen_f(group_id, dataset_name, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be opened."
            stop
        end if

        ! get dataspace
        call h5dget_space_f(dataset_id, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be obtained."
            stop
        end if

        ! read data
        call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, data_shape, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be read."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_read_array_double_3d

    ! The following subroutine reads a 4D array of double precision reals from a HDF5 file
    subroutine hdf5_read_array_double_4d(file_id, group_id, dataset_name, data)
        !> @brief Read a 4D array of double precision reals from a HDF5 file
        !> @param file_id       File id of the file from which the data is read
        !> @param group_id      Group id from which the data is read
        !> @param dataset_name  Name of the dataset from which the data is read
        !> @param data          Data to be read
        !> @author  Ravindra Shinde
        !> @email r.l.shinde@utwente.nl

        integer(hid_t), intent(in)              :: file_id
        character(len=*), intent(in)            :: dataset_name
        double precision, intent(out)           :: data(:,:,:,:)
        integer(hid_t), intent(in)              :: group_id
        integer(hid_t)                          :: dataset_id
        integer(hid_t)                          :: dataspace_id
        integer(HSIZE_T)                        :: data_dims
        integer(HSIZE_T)                        :: data_shape(4)
        integer                                 :: ierr

        data_shape = shape(data)

        ! open dataset
        call h5dopen_f(group_id, dataset_name, dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be opened."
            stop
        end if

        ! get dataspace
        call h5dget_space_f(dataset_id, dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be obtained."
            stop
        end if

        ! read data
        call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, data_shape, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be read."
            stop
        end if

        ! close dataset
        call h5dclose_f(dataset_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataset could not be closed."
            stop
        end if

        ! close dataspace
        call h5sclose_f(dataspace_id, ierr)
        if (ierr /= 0) then
            write(errunit,*) "Error: HDF5 dataspace could not be closed."
            stop
        end if

    end subroutine hdf5_read_array_double_4d

end module hdf5_utils