// Hide the console window on Windows for the GUI build.
#![cfg_attr(
    all(feature = "gui", not(debug_assertions), target_os = "windows"),
    windows_subsystem = "windows"
)]

mod converter;

#[cfg(feature = "gui")]
mod gui;

#[cfg(not(feature = "gui"))]
mod cli;

fn main() {
    #[cfg(feature = "gui")]
    gui::run();

    #[cfg(not(feature = "gui"))]
    cli::run();
}
