def slugify(filename) {
  return filename.toLowerCase().replaceAll(' ', '_').replaceAll('[^a-z0-9._-]', '')
}
