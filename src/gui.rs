use crate::converter::{convert, Options};
use eframe::egui;
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};

struct App {
    src: String,
    dst: String,
    merge: bool,
    scale: String,
    angle: String,
    log: Arc<Mutex<Vec<String>>>,
    busy: Arc<AtomicBool>,
}

impl App {
    fn new(initial: Option<String>) -> Self {
        let mut app = App {
            src: String::new(),
            dst: String::new(),
            merge: true,
            scale: "1.0".to_string(),
            angle: "0.1".to_string(),
            log: Arc::new(Mutex::new(vec![
                "Drop an .stl file onto this window, or click Browse.".to_string(),
                "Flat areas become real planar faces; curved surfaces stay faceted.".to_string(),
            ])),
            busy: Arc::new(AtomicBool::new(false)),
        };
        if let Some(p) = initial {
            app.set_src(PathBuf::from(p));
        }
        app
    }

    fn set_src(&mut self, p: PathBuf) {
        self.dst = p.with_extension("step").to_string_lossy().to_string();
        self.src = p.to_string_lossy().to_string();
    }

    fn push_log(&self, s: &str) {
        self.log.lock().unwrap().push(s.to_string());
    }

    fn start(&mut self, ctx: &egui::Context) {
        if self.src.trim().is_empty() || self.dst.trim().is_empty() {
            self.push_log("ERROR: choose an STL file first.");
            return;
        }
        let scale = match self.scale.trim().parse::<f64>() {
            Ok(v) if v > 0.0 && v.is_finite() => v,
            _ => {
                self.push_log("ERROR: scale must be a positive number.");
                return;
            }
        };
        let angle = match self.angle.trim().parse::<f64>() {
            Ok(v) if v >= 0.0 && v.is_finite() => v,
            _ => {
                self.push_log("ERROR: planar tolerance must be a number >= 0.");
                return;
            }
        };
        let opts = Options {
            merge: self.merge,
            scale,
            angle_tol: angle,
            ..Options::default()
        };
        let src = PathBuf::from(self.src.trim());
        let dst = PathBuf::from(self.dst.trim());
        let log = self.log.clone();
        let busy = self.busy.clone();
        let ctx2 = ctx.clone();
        busy.store(true, Ordering::SeqCst);

        std::thread::spawn(move || {
            let log2 = log.clone();
            let ctx3 = ctx2.clone();
            let mut logger = move |s: String| {
                log2.lock().unwrap().push(s);
                ctx3.request_repaint();
            };
            if let Err(e) = convert(&src, &dst, &opts, &mut logger) {
                logger(format!("ERROR: {e}"));
            }
            busy.store(false, Ordering::SeqCst);
            ctx2.request_repaint();
        });
    }
}

impl eframe::App for App {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // drag & drop an STL onto the window
        let dropped = ctx.input(|i| i.raw.dropped_files.clone());
        if let Some(path) = dropped.iter().find_map(|f| f.path.clone()) {
            self.set_src(path);
        }

        let busy = self.busy.load(Ordering::SeqCst);

        egui::CentralPanel::default().show(ctx, |ui| {
            ui.heading("STL \u{2192} STEP converter");
            ui.add_space(8.0);

            egui::Grid::new("form").num_columns(3).spacing([8.0, 8.0]).show(ui, |ui| {
                ui.label("STL file");
                ui.add(egui::TextEdit::singleline(&mut self.src).desired_width(380.0));
                if ui.button("Browse\u{2026}").clicked() {
                    if let Some(p) = rfd::FileDialog::new()
                        .add_filter("STL files", &["stl", "STL"])
                        .pick_file()
                    {
                        self.set_src(p);
                    }
                }
                ui.end_row();

                ui.label("STEP output");
                ui.add(egui::TextEdit::singleline(&mut self.dst).desired_width(380.0));
                if ui.button("Save as\u{2026}").clicked() {
                    if let Some(p) = rfd::FileDialog::new()
                        .add_filter("STEP files", &["step", "stp"])
                        .save_file()
                    {
                        self.dst = p.to_string_lossy().to_string();
                    }
                }
                ui.end_row();

                ui.label("Scale factor");
                ui.add(egui::TextEdit::singleline(&mut self.scale).desired_width(80.0));
                ui.label("25.4 if the STL is in inches");
                ui.end_row();

                ui.label("Planar tolerance");
                ui.add(egui::TextEdit::singleline(&mut self.angle).desired_width(80.0));
                ui.label("degrees");
                ui.end_row();
            });

            ui.add_space(4.0);
            ui.checkbox(
                &mut self.merge,
                "Merge coplanar triangles into planar faces (recommended)",
            );
            ui.add_space(8.0);

            let width = ui.available_width();
            let btn = ui.add_enabled(
                !busy,
                egui::Button::new(if busy { "Converting\u{2026}" } else { "Convert" })
                    .min_size(egui::vec2(width, 30.0)),
            );
            if btn.clicked() {
                self.start(ctx);
            }

            ui.add_space(8.0);
            let lines: Vec<String> = self.log.lock().unwrap().clone();
            egui::ScrollArea::vertical()
                .stick_to_bottom(true)
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    for l in &lines {
                        ui.monospace(l);
                    }
                });
        });
    }
}

pub fn run() {
    let initial = std::env::args().nth(1);
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([680.0, 540.0]),
        ..Default::default()
    };
    if let Err(e) = eframe::run_native(
        "STL to STEP Converter",
        options,
        Box::new(move |_cc| Ok(Box::new(App::new(initial)))),
    ) {
        eprintln!("GUI error: {e}");
    }
}
