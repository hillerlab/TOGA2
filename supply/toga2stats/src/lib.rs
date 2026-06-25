// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Crate root for `toga2_stats`.
//!
//! CLI argument definitions are in [`cli`] and shared utilities in [`utils`].
//! The binary entry point in `main.rs` ties them together to compute summary
//! statistics from TOGA2 gene-projection output.

pub mod cli;
pub mod utils;
