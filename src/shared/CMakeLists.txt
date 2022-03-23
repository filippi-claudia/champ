# -----------------------------
# These are the objects files which are both present in dmc and vmc
# Depends on shared_modules
# -----------------------------

file(GLOB SHARED_OBJECTS_F90_ "*.f90")
file(GLOB SHARED_OBJECTS_F_ "*.f")

# --------------------------------------------------------------------
# Remove objects to share depending on CMAKE options:
# --------------------------------------------------------------------
if (ENABLE_QMMM)
  list(REMOVE_ITEM SHARED_OBJECTS_F_
    ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_placeholder.f)
else()
  list(REMOVE_ITEM SHARED_OBJECTS_F90_
    ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_bspline.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_modules.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_pot.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_density.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_writecube.f90)
endif()

if (ENABLE_PERIODIC)
  list(REMOVE_ITEM SHARED_OBJECTS_F_
    ${CMAKE_CURRENT_SOURCE_DIR}/pw_placeholder.f)
else()
  list(REMOVE_ITEM SHARED_OBJECTS_F_
    ${CMAKE_CURRENT_SOURCE_DIR}/pw_ewald.f
    ${CMAKE_CURRENT_SOURCE_DIR}/pw_orbitals_e.f
    ${CMAKE_CURRENT_SOURCE_DIR}/pw_find_image.f
    ${CMAKE_CURRENT_SOURCE_DIR}/pw_orbitals.f
    ${CMAKE_CURRENT_SOURCE_DIR}/pw_read.f)
endif()

set(SHARED_OBJECTS_F ${SHARED_OBJECTS_F_} PARENT_SCOPE)
set(SHARED_OBJECTS_F90 ${SHARED_OBJECTS_F90} PARENT_SCOPE)
set(SHARED_OBJECTS ${SHARED_OBJECTS_F90_} ${SHARED_OBJECTS_F_} PARENT_SCOPE)

