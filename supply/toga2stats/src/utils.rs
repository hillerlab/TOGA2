// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Shared utility functions.
//!
//! Currently provides [`reader`], a simple convenience wrapper around
//! synchronous file I/O.

use std::fmt::Debug;
use std::fs::File;
use std::io::Read;
use std::path::Path;

/// Reads the entire contents of a file into a `String`.
///
/// Generic over any path-like type that implements `AsRef<Path>` and `Debug`.
///
/// # Arguments
///
/// * `file` - Path to the file to read.
///
/// # Returns
///
/// A `Result<String, Box<dyn std::error::Error>>` containing the file contents
/// on success, or an error if the file cannot be opened or read.
pub fn reader<P: AsRef<Path> + Debug>(file: P) -> Result<String, Box<dyn std::error::Error>> {
    let mut file = File::open(file)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    Ok(contents)
}
