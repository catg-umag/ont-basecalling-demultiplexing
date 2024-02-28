/*
 * Check if file exists of throw error
 */
def pathCheck(path, isDirectory = false) {
  f = file(path)
  if (!f.exists()) {
    throwError("Path ${path} does not exist")
  } else if (!f.isFile() && !isDirectory) {
    throwError("Path ${path} is not a file")
  } else if (!f.isDirectory() && isDirectory) {
    throwError("Path ${path} is not a directory")
  }

  return f
}

def sanitizeFilename(filename) {
  return filename.toLowerCase().replaceAll(' ', '_').replaceAll('[^a-z0-9._-]', '')
}

def throwError(msg) {
  exit 1, "ERROR: $msg"
}