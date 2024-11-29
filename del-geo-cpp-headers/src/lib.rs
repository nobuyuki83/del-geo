pub struct Headers {}

pub const HEADERS: Headers = Headers {};

impl Headers {
    pub fn get(&self, idx: u32) -> Option<(&str, &str)> {
        match idx {
            0 => Some(("aabb.h", include_str!("aabb.h"))),
            1 => Some(("aabb2.h", include_str!("aabb2.h"))),
            2 => Some(("aabb3.h", include_str!("aabb3.h"))),
            3 => Some(("vec3.h", include_str!("vec3.h"))),
            4 => Some(("tri3.h", include_str!("tri3.h"))),
            5 => Some(("mat2_sym.h", include_str!("mat2_sym.h"))),
            6 => Some(("mat2x3_col_major.h", include_str!("mat2x3_col_major.h"))),
            7 => Some(("mat3_col_major.h", include_str!("mat3_col_major.h"))),
            8 => Some(("mat4_col_major.h", include_str!("mat4_col_major.h"))),
            _ => None,
        }
    }

    pub fn write_files(&self, path_out_dir: &std::path::PathBuf) {
        for header_idx in 0.. {
            use std::io::Write;
            let Some((header_name, header_content)) = self.get(header_idx) else {
                break;
            };
            let path = path_out_dir.join(header_name.to_string());
            let mut output = std::fs::File::create(path).unwrap();
            write!(output, "{}", header_content).unwrap();
        }
    }
}
