import math
import itk
import matplotlib.pyplot as plt

input_name = "D:\\skola\\DP\\dp_isr\\sample_files\\dicom\\"

# rotation angles of transformation (degrees)
rx = -90.0
ry = 0.0
rz = 0.0

# translation values of transformation (mm)
tx = 0.0
ty = 0.0
tz = 0.0

# center of transformation (mm)
cx = 0.0
cy = 0.0
cz = 0.0

# distances from DICOM metadata
sid = 1085.6  # source to detector distance
spd = 595.0  # source to patient distance
# sx = sy = 1.0

# output DRR image sizes
dx = 512.0
dy = 779.0

# output DRR image origin location
# o2Dx = 0.0
# o2Dy = 0.0

# Hounsfield values below this threshold not taken into
# account during X-ray simulation
# threshold of 200 reasonably omits soft tissue
threshold = 200

dimension = 3

pixel_type = itk.UC
input_pixel_type = itk.SS
output_pixel_type = itk.UC
input_image_type = itk.Image[input_pixel_type, dimension]
output_image_type = itk.Image[output_pixel_type, dimension]

# getting DICOM files and loading CT volume
image_type = itk.Image[input_pixel_type, dimension]
reader_type = itk.GDCMSeriesFileNames.New()
reader_type.SetUseSeriesDetails(True)
reader_type.SetDirectory(input_name)
file_names = reader_type.GetInputFileNames()

reader = itk.ImageSeriesReader[image_type].New()
dicom_io = itk.GDCMImageIO.New()
dicom_io.LoadPrivateTagsOn()
reader.SetImageIO(dicom_io)
reader.SetFileNames(file_names)
reader.Update()
image = reader.GetOutput()
metadata = dicom_io.GetMetaDataDictionary()

print(f"region: {image.GetBufferedRegion()}, resolution: {image.GetSpacing()}, origin: {image.GetOrigin()}, "
      f"type: {type(image)}")

# getting center of input volume to rotate x-ray - DRR image ("detector") around
im_center = image.GetBufferedRegion().GetSize()
im_center[0] //= 2
im_center[1] //= 2
im_center[2] //= 2 # index values in pixels
im_center_mm = image.TransformIndexToPhysicalPoint(im_center) # output in milimeters

# filter for the output DRR image
resample_filter = itk.ResampleImageFilter[input_image_type, input_image_type].New()
resample_filter.SetInput(image)
resample_filter.SetDefaultPixelValue(0)

# transform for the x-ray source and DRR image ("detector")
transform = itk.CenteredEuler3DTransform[itk.D].New()
transform.SetComputeZYX(True)
translation = [tx, ty, tz]

deg_to_rad = (math.atan(1.0) * 4.0) / 180.0
transform.SetRotation(deg_to_rad * rx, deg_to_rad * ry, deg_to_rad * rz)
im_origin = image.GetOrigin()
im_resolution = image.GetSpacing()
im_region = image.GetBufferedRegion()
im_size = im_region.GetSize()

# im_origin[0] += im_resolution[0] * im_size[0] / 2.0
# im_origin[1] += im_resolution[1] * im_size[1] / 2.0
# im_origin[2] += im_resolution[2] * im_size[2] / 2.0
center = [cx + im_center[0], cy + im_center[1], cz + im_center[2]]
transform.SetCenter(center)
print(f"Image Size: {im_size}\n"
      f"Resolution: {im_resolution}\n"
      f"origin: {im_origin}\n"
      f"center: {center}\n"
      f"Transform: {transform}\n")

# creating x-ray source and setting its location - focal point
interpolator = itk.itkRayCastInterpolateImageFunctionPython.itkRayCastInterpolateImageFunctionISS3D.New()
interpolator.SetTransform(transform)
interpolator.SetThreshold(threshold)

# focal point's location in x and y-axis same as volume origin,
# moved away from input volume along z-axis by source to patient distance
focal_point = [im_origin[0],
               im_origin[1],
               im_origin[2] - spd]
interpolator.SetFocalPoint(focal_point)

print(f"Focal point: {focal_point}\n"
      f"Interpolator: {interpolator}")

resample_filter.SetInterpolator(interpolator)
resample_filter.SetTransform(transform)
size = itk.Size[dimension]((int(dx), int(dy), 1))
resample_filter.SetSize(size)

# output DRR image spacing same as volume's spacing and slice thickness
spacing = [im_resolution[0], im_resolution[2], 1.0]
resample_filter.SetOutputSpacing(spacing)

print(f"Output image size: {size}\n"
      f"Output image spacing: {spacing}")

# setting output DRR image location
# z-axis: distance between x-ray source and DRR image is source to detector distance
origin = [focal_point[0],
          focal_point[1],
          focal_point[2] + sid]

resample_filter.SetOutputOrigin(origin)

print(f"Output image origin: {origin}")

# rescaling output intensity values and casting to uint8
rescaler = itk.RescaleIntensityImageFilter[input_image_type, input_image_type].New()
rescaler.SetInput(resample_filter.GetOutput())
rescaler.SetOutputMinimum(0)
rescaler.SetOutputMaximum(255)

cast_image_filter = itk.CastImageFilter[input_image_type, output_image_type].New()
cast_image_filter.SetInput(rescaler.GetOutput())
cast_image_filter.Update()

# displaying the result
image_view = itk.GetArrayViewFromImage(cast_image_filter.GetOutput())

plt.imshow(image_view[0], cmap="gray")
plt.axis("off")
plt.show()
