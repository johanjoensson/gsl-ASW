#ifndef SITE_H
#define SITE_H
#include "GSLpp/vector.h"

template<size_t dim>
class Site_t;

namespace std {
	template<size_t dim>
	struct hash<Site_t<dim>>{
		hash() = default;
		size_t operator()(Site_t<dim> &s) const
		{
			size_t res = 0;
			for(const auto& elem : s.pos()){
				res ^= std::hash<double>()(elem);
			}
			return (res << 1) ^ (std::hash<size_t>()(s.index()) >> 1);
		}
	};
}

template<size_t dim>
class Site_t{
	private:
		size_t index_m;
		std::array<size_t, dim> coord_m;
		GSL::Vector pos_m;

		Site_t(size_t idx, std::array<size_t, dim> c, GSL::Vector&& p)
		 : index_m(std::move(idx)), coord_m(std::move(c)), pos_m(std::move(p))
		{}
	public:
		// Site_t() : index_m(), coord_m(), pos_m() {}
		Site_t(const size_t idx, const GSL::Vector& pos, const std::array<size_t, dim>& size)
		 : index_m(idx), coord_m(calc_coord(idx, size)), pos_m(pos.clone())
		{}

		Site_t(const std::array<size_t, dim>& coords, const GSL::Vector& pos, const std::array<size_t, dim>& size)
		 : index_m(calc_index(coords, size)), coord_m(coords), pos_m(pos.clone())
		{}


		Site_t clone() const
		{
			return Site_t(this->index_m, this->coord_m, this->pos_m.clone());
		}

		std::array<size_t, dim> calc_coord(const size_t index, const std::array<size_t, dim>& size)
		{
			size_t tmp = index;
			std::array<size_t, dim> res;
			for(size_t i = dim; i > 0; i--)
			{
				size_t offset = 1;
				for(size_t j = 0; j < i - 1; j++ ){
					offset *= size[j];
				}
				res[i - 1] = tmp / offset;
				tmp -= res[i - 1]*offset;
			}
			return res;
		}

		size_t calc_index(const std::array<size_t, dim>& coords, const std::array<size_t, dim>& size)
		{
			size_t res = 0;
			for(size_t i = dim; i > 0; i--)
			{
				size_t offset = 1;
				for(size_t j = 0; j < i - 1; j++ ){
					offset *= size[j];
				}
				res += offset * coords[i - 1];
			}
			return res;
		}

		size_t index() const {return index_m;}
		const std::array<size_t, dim>& coord() const {return coord_m;}
		GSL::Vector::Const_View pos() const {return pos_m.cview();}
		// const GSL::Vector pos() const {return pos_m.view();}
		void set_pos(const GSL::Vector& pos){pos_m.copy(pos);}

		bool operator==(const Site_t<dim>& s) const
		{
			return (index() == s.index() && coord() == s.coord() && pos() == s.pos());
		}

		bool operator!=(const Site_t<dim>& s) const
		{
			return !(*this == s);
		}
};

template<size_t dim>
struct Site_t_hasher
{
	size_t operator()(const Site_t<dim>& site) const
	{
		size_t res = 0;
		for(const auto& elem : site.coord()){
			res ^= std::hash<double>()(elem);
		}
		res <<= 1;
		for(const auto& elem : site.pos()){
			res ^= std::hash<double>()(elem);
		}
		return std::hash<size_t>()(site.index())^(res >> 1);
	}
};
#endif // SITE_H
