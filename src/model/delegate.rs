macro_rules! delegate_analytical_material {
    ($type:ident) => {
        impl<R> lamina_core::Material for $type<R>
        where
            R: num_traits::Float + std::fmt::Debug,
        {
            type Real = R;

            fn relative_permeability<I, C>(&self, vacuum_wavenumber: I) -> I::Mapped<C>
            where
                C: lamina_core::ComplexScalar<RealField = R> + Copy,
                I: lamina_core::Sampled<Elem = R>,
            {
                self.inner.relative_permeability(vacuum_wavenumber)
            }

            fn relative_permittivity<I, C>(&self, vacuum_wavenumber: I) -> I::Mapped<C>
            where
                C: lamina_core::ComplexScalar<RealField = R> + Copy,
                I: lamina_core::Sampled<Elem = R>,
            {
                self.inner.relative_permittivity(vacuum_wavenumber)
            }
        }

        impl<R> lamina_core::DifferentiableMaterial for $type<R>
        where
            R: num_traits::Float + std::fmt::Debug,
        {
            fn relative_permittivity_derivative<I, C>(
                &self,
                vacuum_wavenumber: I,
                order: lamina_core::DerivativeOrder,
            ) -> I::Mapped<C>
            where
                C: lamina_core::ComplexScalar<RealField = R> + Copy,
                I: lamina_core::Sampled<Elem = R>,
            {
                self.inner
                    .relative_permittivity_derivative(vacuum_wavenumber, order)
            }
        }

        impl<R> lamina_core::MeromorphicMaterial for $type<R>
        where
            R: num_traits::Float + std::fmt::Debug,
        {
            fn relative_permeability_complex<I, C>(&self, vacuum_wavenumber: I) -> I::Mapped<C>
            where
                C: lamina_core::ComplexScalar<RealField = R> + Copy,
                I: lamina_core::Sampled<Elem = C>,
            {
                self.inner.relative_permeability_complex(vacuum_wavenumber)
            }

            fn relative_permittivity_complex<I, C>(&self, vacuum_wavenumber: I) -> I::Mapped<C>
            where
                C: lamina_core::ComplexScalar<RealField = R> + Copy,
                I: lamina_core::Sampled<Elem = C>,
            {
                self.inner.relative_permittivity_complex(vacuum_wavenumber)
            }
        }

        impl<R> lamina_core::DifferentiableMeromorphicMaterial for $type<R>
        where
            R: num_traits::Float + std::fmt::Debug,
        {
            fn relative_permittivity_complex_derivative<I, C>(
                &self,
                vacuum_wavenumber: I,
                order: lamina_core::DerivativeOrder,
            ) -> I::Mapped<C>
            where
                C: lamina_core::ComplexScalar<RealField = R> + Copy,
                I: lamina_core::Sampled<Elem = C>,
            {
                self.inner
                    .relative_permittivity_complex_derivative(vacuum_wavenumber, order)
            }
        }
    };
}

pub(crate) use delegate_analytical_material;
