pub struct Headers {}

pub const HEADERS: Headers = Headers {};

impl Headers {
    pub fn get(&self, idx: usize) -> Option<(&str, &str)> {
        [
            ("aabb.h", include_str!("aabb.h")),
            ("aabb2.h", include_str!("aabb2.h")),
            ("aabb3.h", include_str!("aabb3.h")),
            ("mat2_sym.h", include_str!("mat2_sym.h")),
            ("mat2x3_col_major.h", include_str!("mat2x3_col_major.h")),
            ("mat3_col_major.h", include_str!("mat3_col_major.h")),
            ("mat4_col_major.h", include_str!("mat4_col_major.h")),
            ("quaternion.h", include_str!("quaternion.h")),
            ("tri3.h", include_str!("tri3.h")),
            ("vec3.h", include_str!("vec3.h")),
            ("edge2.h", include_str!("edge2.h")),
            ("tri2.h", include_str!("tri2.h")),
            ("vec2.h", include_str!("vec2.h")),
        ]
        .get(idx)
        .copied()
    }

    pub fn write_files(&self, path_out_dir: &std::path::Path) {
        for header_idx in 0.. {
            use std::io::Write;
            let Some((header_name, header_content)) = self.get(header_idx) else {
                break;
            };
            let path = path_out_dir.join(header_name);
            let mut output = std::fs::File::create(path).unwrap();
            write!(output, "{}", header_content).unwrap();
        }
    }
}
